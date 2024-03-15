# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2022 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    Module commandprocessor.  Main handling of subprocess calls -- currently all are bash shell runs.
    Will eventually expand to async multiprocessor handling with a list of processes.  To allow queue'ed calls and
    the main program window to be non-blocking.  Handles new ProcessWait pop-up with timer to terminate errant
    subprocess calls. Still not terminate for running call though.

    v1 /v2 was 50% BATCH and 50% BASH calls.  Moving to all BASH, code is more OS
    independent (as long as Win10 environment has BASH/Unix calls available) and internal file / path handling easier
    as most are left in universal/Unix format and not nativeOS specific. Only the BASH script filename itself needs
    to be in nativeOS form. Historically, BATCH files were needed as dos2unix and similar calls were made.  But all
    dos2unix needs were replaced by simply writing binary ('wb') and using .encode() on strings written in binary
    form; alleviating the need to convert back to common POSIX \n format.
"""

import os           # os.chmod()
import stat         # stat.S_IXxxx flags
import time         # time.time()
import subprocess   # Popen, run, etc
import platform
if platform.uname().system != "Linux":
    from wakepy import set_keepawake, unset_keepawake   # Ugly, but cannot import on Ubuntu aparently
else:
    # from elevate import elevate                       # Only needed for wakepy on Linux
    # from wakepy import set_keepawake, unset_keepawake
    pass

from tkinter import Toplevel, Label

from utilities import DEBUG, wgse_message
import settings as wgse
# from mainwindow import mainwindow_resume     # embedded in button_continueAfterBatchJob to break loop


####################################################################################################################
# Command Execution Subsection

pleaseWaitWindow = None


def is_command_available(command, opt, internal=True):
    """ Verify a shell command is available to run as specified.  Get version to check as well. """
    import re

    # Table of known program version return strings that we can extract a version from; expect in first line returned
    # --version on 11+ gives non-quoted, without "version". -version on all gives as shown except 11+ start with 11+
    lookup = {
      'openjdk v': r"openjdk version \"1[.](?P<major>\d+)[.](?P<minor>\d+)_(?P<sub>\d+)",  # openjdk version "1.8.0_312"
      'openjdk 1': r"openjdk (?P<major>\d+)[.](?P<minor>\d+)[.](?P<sub>\d+)",              # openjdk 11.0.14 2022-02-08
      'openjdk 2': r"openjdk (?P<major>\d+)[.](?P<minor>\d+)[.](?P<sub>\d+)",              # openjdk 20.x.xx 2023-03-21
      'java vers': r"java version \"1[.](?P<major>\d+)[.](?P<minor>\d+)_(?P<sub>\d+)",     # Oracle Java 8
      'samtools ': r"samtools (?P<major>\d+)[.](?P<minor>\d+)[.]?(?P<sub>[a-z0-9]*)",      # samtools 1.12  or 1.15.1
      'python 3.': r"python (?P<major>\d+)[.](?P<minor>\d+)[.](?P<sub>\d+)"                # python 3.8.9
    }

    try:
        result = subprocess.run([command, opt], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except:
        if not internal:
            wgse_message("warning", 'NoCmdWindowTitle', True,
                         wgse.lang.i18n['NoCmdInstalledErrorMessage'].replace("{{cmd}}", command))
        return False
    else:
        # Ugly; when True we pass back the major / minor version as a global in settings.py
        cp_version = ['0', '0', '0']    # Start with list of strings; end with list of ints
        check = result.stdout.decode('utf-8', errors='ignore').splitlines()     # Windows WSL2 generates errors
        if check and len(check) and len(check[0]) > 9:
            # Use the 1st 9 characters of 1st command return line to lookup regex string to extract major/minor version
            check[0] = check[0].lower()
            pattern = re.compile(lookup.get(check[0][:9], '(?P<major>garbage).(?P<minor>trash).(?P<sub>recycle)'))
            match = pattern.match(check[0])         # Find version regex string match
            if match:
                cp_version = match.groups("0")
        else:
            if not internal:
                wgse_message("warning", 'NoVersionTitle', True,
                             wgse.lang.i18n['NoVersionMessage'].replace("{{cmd}}", command))
        if len(cp_version[2]) == 0:     # Special case when only two elements
            cp_version = [cp_version[0], cp_version[1], '0']
        wgse.cp_version = list(map(int, cp_version))
        return True


def simple_command(command_and_opts):
    """
    Intended for simple shell-level commands that return a single value / line.  Simply runs the command and returns
    the result.  For especially when the value cannot be obtained from a python library call. sysctl, cpuinfo, etc.
    Command_and_opts must be list of strings of each opt for the first string which is the command.
    """
    try:
        result = subprocess.run(command_and_opts, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except:
        return None
    else:
        return result.stdout.decode('utf-8', errors='ignore').strip()   # Windows WSL2 generates errors


def time_label(etime):
    time_exp = (etime / 3600., 'hours')      if etime > 3600 else \
                (int(etime / 60), 'minutes') if etime > 60 else \
                (etime, 'seconds')
    return f"{time_exp[0]:3.1f} {time_exp[1]}" if time_exp[1] == 'hours' else f"{time_exp[0]:3d} {time_exp[1]}"


def run_external_program(script_FBS, command):
    """ Run an external batch program with a time-limit"""

    # Due to adding the "Direct" mode to run_bash_script, script may not be last entry in command list
    command_str = os.path.basename(command[-1] if ".sh" in command[-1] else command[1])  # " ".join(command).strip()

    # Rough std deviation (max) from expected time (x2); more robust for table out-of-date
    maxtime = wgse.expected_time.get(script_FBS, 3600) * 2

    start_time = time.time()
    start_ctime = time.ctime()

    print(f'--- Exec: {command_str}, started @ {start_ctime}')     # REH 14Mar2020 Debug
    # Blocking subprocess call (by using p.comminucate() )
    try:
        # subprocess.run(command_to_run_and_args, timeout=maxtime*2)
        with subprocess.Popen(command) as p:  # , stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True) as p:
            # stdout, stderr = p.communicate(timeout=maxtime)
            stdout, stderr = p.communicate()
    except subprocess.TimeoutExpired:   # Todo This specific exception not being caught here?
        # Todo need to see if due to a timeout and handle appropriately; also user initiated abort
        p.kill()
        stdout, stderr = p.communicate()
        stop_ctime = time.ctime()
        tlabel = time_label(maxtime)
        DEBUG(f"--- FAILURE to finish before timeout {tlabel}: {command_str} (@ {stop_ctime})----")
        raise
    except subprocess.CalledProcessError:
        # Todo need to handle different process termination conditions
        DEBUG(f"--- FAILURE to execute: {command_str} ----")
        raise
    else:
        stop_ctime = time.ctime()
        tlabel = time_label(round(time.time() - start_time))
        print(f'--- SUCCESS: {tlabel} to run: {command_str} (finished @ {stop_ctime}')


def show_please_wait_window(reason, parent):
    global pleaseWaitWindow

    # Make code more robust to errors; external language file and expected time table may not be up-to-date
    etime = wgse.expected_time.get(reason, 0)
    # etime = round(min(1.5, etime * wgse.BAM.relfsize)) # See note below on adjusting time relative to 60 GB nominal
    tlabel = time_label(etime) if etime else f"unknown ({reason}?)"
    verbose_reason = wgse.lang.i18n.get(reason, f"(Internal Translation Error: {reason} unknown)")
    DEBUG(f'In Please Wait: {verbose_reason}, Expected Wait: {tlabel}, Start')
    """
        We will adjust expected wait time by amount relative to nominal 45 GB 30x WGS file size. Floor of 0.01 adj.
        Using linear but likely not linear.  Idea is smaller, subsetted BAMs may take less time? True for operations
        that have to scan the whole file but not others?  Need to study this further and drop if nonsense.
        Does not seem to work as written; commenting out for now
    """

    # Due to a limitation in the keepawake on Linux, it requires elevated (SUDO) privileges. we can request the sudo
    # and fork the new job in elevated mode. But ugly as it creates a new command window. And stays on until WGSE is
    # exited. No way to request a drop elevated. So just ignore and document that the feature does not work on Linux.
    if wgse.os_plat == "Linux":     # and etime > 3600 and os.getuid() != 0:
        # wgse.DEBUG_MODE = True  # We turn on Debug Mode so some of the open-ended "rm -rf" are not executed under sudo
        # elevate(show_console=True)  # Request elevated privilages from the user
        # set_keepawake()
        pass
    else:
        set_keepawake()

    font = wgse.fonts.table if wgse and wgse.fonts else \
        {'14': ("Times New Roman", 14), '28b': ("Arial Black", 28, "Bold")}

    pleaseWaitWindow = Toplevel(wgse.window)
    pleaseWaitWindow.transient(parent)
    pleaseWaitWindow.title(wgse.lang.i18n['PleaseWait'])
    pleaseWaitWindow.geometry("")
    pleaseWaitWindow.protocol("WM_DELETE_WINDOW", cancelWait)
    # pleaseWaitWindow.attributes('-topmost', 'true')
    # pleaseWaitWindow.parent.withdraw()
    pleaseWaitWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    pleaseWaitWindow.columnconfigure(0, weight=1)
    pleaseWaitWindow.rowconfigure(0, weight=1)
    Label(pleaseWaitWindow, text=str(wgse.lang.i18n['PleaseWait']), font=font['28b']).grid(column=0, row=0, padx=56, pady=28)
    Label(pleaseWaitWindow, text=verbose_reason, font=font['14']).grid(column=0, row=1, padx=1, pady=1)
    Label(pleaseWaitWindow, text=f'{wgse.lang.i18n["ExpectedWait"]} {tlabel}.').grid(column=0, row=2, padx=1, pady=1)
    time_now = time.ctime()
    Label(pleaseWaitWindow, text=f'{wgse.lang.i18n["StartedAt"]} {time_now}').grid(column=0, row=3, padx=1, pady=1)
    # TODO add user abort button that is caught, kills job waiting on, and takes one back where?
    pleaseWaitWindow.update()
    pleaseWaitWindow.grab_set()


def cancelWait():
    global pleaseWaitWindow

    # DEBUG("cancelWait; unfortunately not yet implemented!")
    # TODO Setup exception that can be caught by threads executing pleaseWait
    # noinspection PyUnresolvedReferences
    if wgse.os_plat != "Linux":
        unset_keepawake()
    pleaseWaitWindow.destroy()


def abortWait():
    # Todo implement abortWait()
    # DEBUG("abortWait called; not yet implemented.")
    cancelWait()


def finishWait():
    # Todo implement finishWait()
    # DEBUG("finishWait called; not yet implemented.")
    cancelWait()


def run_bash_script(script_title, script_contents, parent=wgse.window, direct=False):
    """
        Main entry point for commandprocessor module.  Two modes though: Direct or not.  In Direct, the
        script_contents is a shlex like list of the single, parsed command line.  No need for quotes, etc.
        In not Direct, we create a BASH script file in Temp from the supplied (multi-line) command content
        in a string, then run the newly created bash script directly.
    """
    if not direct:
        script_oFN = f'{wgse.tempf.oFP}{script_title}.sh'

        if os.path.isfile(script_oFN):
            os.remove(script_oFN)

        with open(script_oFN, "wb") as f:          # 15 Mar 2020 REH Changed to binary write to keep Unix \n format
            f.write("#!/usr/bin/env bash -x\n".encode())
            f.write(script_contents.encode())
        os.chmod(script_oFN, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # Historic; needed as sourced by BASH now?

        # Full path specified so do not need to pre-pend './' for Linux; Windows needs BASH command to start .sh files
        # BASH script file. So no need to parse script content using shlex. We create our own parsed list here
        command = [wgse.bashx_oFN, "-x", script_oFN]  # if wgse.os_plat == "Windows" else [script_oFN]
    else:
        command = script_contents   # Already passed as a list; not a multi-line string like above

    DEBUG(f'Starting command: {" ".join(command).strip()}')

    if wgse.window and wgse.dnaImage:       # if windowing system created and setup (for outside, command line calls)
        show_please_wait_window(script_title, parent)

    run_external_program(script_title, command)

    if wgse.window and wgse.dnaImage:
        finishWait()
