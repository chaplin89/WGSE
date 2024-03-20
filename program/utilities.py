# coding: utf8
#
# Main utilities support
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
  General Utilities module support for the WGS Extract system
    This file is included by all others. To avoid loops, some imports are buried in the functions where needed

  Classes: TemporaryFiles, LanguageStrings, OutputDirectory, FontTypes
  Routines to manipulate JSON files, and upgrade the WGSE program itself based on JSON version files
"""
import os           # for os.path.*, os.remove, os.listdir
import shutil       # shutil.rmtree(), shutil.diskusage

import settings as wgse


class Error(Exception):
    pass


class Warning(Exception):
    pass


def rreplace(pat, sub, string):
    return string[:-len(pat)] + sub if string.endswith(pat) else string


def bam_base(sam_file):
    return rreplace(".bam", "", rreplace(".cram", "", rreplace(".sam", "", sam_file)))


def fastq_base(fastq_file):
    # Remove possible .gz before regular extension
    return rreplace(".fastq", "", rreplace(".fq", "", rreplace(".gz", "", fastq_file)))


def fasta_base(fasta_file):
    # Remove possible .gz before regular extension
    return rreplace(".fasta", "", rreplace(".fa", "", rreplace(".fna", "", rreplace(".gz", "", fasta_file))))


def vcf_base(vcf_file):
    # Remove possible .gz before regular extension
    return rreplace(".vcf", "", rreplace(".bcf", "", rreplace(".gz", "", vcf_file)))


def time_label(etime):
    time_exp = (etime / 3600., 'hours')      if etime > 3600 else \
                (int(etime / 60), 'minutes') if etime > 60 else \
                (etime, 'seconds')
    return f"{time_exp[0]:3.1f} {time_exp[1]}" if time_exp[1] == 'hours' else f"{time_exp[0]:3d} {time_exp[1]}"


def DEBUG(msg):
    """ Old C Macro hack for programmer debug system; but remains as executed code in Python """
    if wgse and wgse.DEBUG_MODE:         # Could use __debug__ but does not save much as calls have if statement then
        print(f"DEBUG: {msg}")


def offset(root_geometry, amount):
    """ For use in tkinter sub-windows; to add an offset to the root windows location for the geometry command.
        Takes the result of the root.geometry() command that includes size and then position. Returns just
        the position string that can be appended to a desired size string or used as is.
    """
    return "+" + "+".join([str(int(i) + amount) for i in root_geometry.split('+')[1:]])


def wgse_message(etype, title, bodyx, body):
    """ Global message utility.  But only show in interactive / GUI mode or if DEBUG is on"""
    from mainwindow import gui_popup

    if wgse.lang and wgse.lang.i18n:
        ttitle = wgse.lang.i18n[title]
        tbody = body if bodyx else wgse.lang.i18n[body]    # Is body translated yet? If not, then translate
        ptbody = tbody.replace("\n\n", "\n")
    else:
        ttitle = title
        tbody = body + ' (no lang set)'
        ptbody = body

    if wgse.gui or wgse.DEBUG_MODE:
        if etype == "error":
            print(f'***ERROR: {ptbody}')
        elif etype == "warning":
            print(f'---WARN: {ptbody}')
        elif etype == "info":
            print(f'+++INFO: {ptbody}')
        # Other etypes (yesno, okcancel, yesnocancel) are all questions and only handled by the gui popup in the mainwindow

    if wgse.gui and wgse.window:
        return gui_popup(etype, ttitle, tbody)      # Some are questions with return of None, True and/or False



# Technically, Cygwin BASH will accept both forms.  And in fact, WinPython returns Win10 disk drive letters to
# start a path. Which most programs (CygWin bash included) seem to still accept even if they require a forward slash.
def nativeOS(path_FP):
    """
    From Universal to native OS specific.  Replace os_slash forward slash to back slash if Win10.
    Also handles drive designations buried in pseudo-Unix paths to proper Win10 system paths.
    Nothing to do on non-Win10 systems.
    """
    if wgse.os_plat == "Windows" and path_FP:
        if path_FP[0:9] == "/cygdrive/":            # Just in case; should never occur
            path_FP = f'{path_FP[10].upper()}:{path_FP[11:]}'    # Convert to Win10 Drive letter format
        elif path_FP[0:4] == "/mnt/":               # In case used WSL mechanism
            path_FP = f'{path_FP[5].upper()}:{path_FP[6:]}'      # Convert to Win10 Drive letter format
        elif path_FP[1] == ":":                    # Often drive letters still come through on Win10 systems
            pass
        elif path_FP[0:1] == "//":                 # Network drive
            pass
        path_oFP = path_FP.replace("/", "\\")       # Change to backward slash
    else:
        path_oFP = path_FP      # Nothing to do for Linux, Unix, MacOS
    # DEBUG(f'nativeOS: before {path_FP}, after {path_oFP}')
    return path_oFP


def universalOS(path_oFP, wsl_fix=False):
    """
    From Native to Universal OS Specific (back slash to forward slash on Win10)
    Is a no-op except in Win10.  Win10 has for forward / back slash change.
    wsl bwa is the only command line that will NOT accept win10 drive letters; so remove drive letters here
    """
    if wgse.os_plat == "Windows" and path_oFP:
        path_FP = path_oFP.replace("\\", "/")      # Change to forward slash
        if path_FP[1] == ":":                      # Win10 Drive letter
            if wsl_fix:
                path_FP = f'/mnt/{path_FP[0].lower()}{path_FP[2:]}'  # Change driver letter to path
            else:
                # path_FP = f'/cygdrive/{path_FP[0].lower()}{path_FP[2:]}'   # Change driver letter to path
                # Leave as is; most programs do not follow or understand Cygwin introduced mount points
                # and Cygwin Bash understands windows disk letter specifiers
                pass
        elif path_FP[0:1] == "//":                  # Network drive; leave alone
            pass
    else:
        path_FP = path_oFP                          # Nothing to do for Linux, Unix, MacOS
    # DEBUG(f'universalOS: before {path_oFP}, after {path_FP}')
    return path_FP


def unquote(path):
    """ Simple routine to remove quotes surrounding a path; instead of simply calling .replace('"', '') """
    return path.replace('"', '') if path else None


# Todo Still needed?  Historical, before quoted path names everywhere.  Added space as legal character now.
def is_legal_path(path_to_check):
    """
    Check file path string for illegal characters. Space no longer illegal as have quoted paths everywhere.
    (but with quoted paths is this check even still needed? Or needed outside of program startup check?)
    """
    import re
    return re.match(r'^[\w%/\-:~.\s\\]+$', path_to_check)


def check_exists(file_oFN, error_mesg):
    import os       # Utilities gets imported all over so hide this import inside routine

    exists = os.path.exists(file_oFN)
    if not exists:
        wgse_message('error', 'errNoFileExists', True, wgse.lang.i18n[error_mesg].replace("{{FILE}}", file_oFN))
    return exists


######################################################################################################################
# Todo Temp Location should become part of a settings class; that also saves / restores settings
class TemporaryFiles:
    """
    A class to manage temporary files.
    Includes self.list for list of files and folders to clean; likely in the user specified output directory
    Optionally clears out user specified temp directory of all files.

    As of 1 Apr 2021, we have modified to use the os.getpid() of the process instance to create a subdirectory
    within the specified Temporary Files directory.  That way, multiple instances of the program can be safely
    started and run.
    """

    def __init__(self, default_FP, top_level=False):

        self.set  = None
        self.FP   = None
        self.oFP  = None
        self.FB   = None            # Now will always be the PID directory created inside tempf
        self.free = 0               # Disk free space in temporary directory area

        self.default_FP = None      # For the root "default"
        self.save_FP = None

        self.change(default_FP)     # Initial call uses this value to define default_FP

        self.list = []  # files to be cleaned that are not inside TempFiles area (helps limit os.remove calls)
        if self.FP:
            self.clean(True) if top_level else ""  # Start with a clean temporary files directory
        else:
            # Setting temporary files directory failed; really fatal and should probably exit (raise exception)
            DEBUG(f'***ERROR: Cannot set default Temporary Files directory: {default_FP}')
            del self

    def change(self, dir_FP):       # For when changing temporary files directory; including initialization
        """
        The real guts of the Temporary Files Directory setup.

        Note: no longer require the directory be empty.  We will setup a subdirectory based on the PID which will,
        by definition, be empty when we start.  OK if other files in the directory to start. This way we now support
        multiple instances of the program at the same time.
        """

        if not dir_FP:
            # Parameter with new Temporary Dirctory must exist
            return      # Simply return and ignore as null may be due to bad restore settings or user cancel

        dir_FP += '/' if dir_FP[-1] != '/' else ''  # Assure trailing slash; universalOS so always forward slash
        dir_oFP = nativeOS(dir_FP)

        if not (is_legal_path(dir_FP) and os.path.isdir(dir_oFP)) or \
           dir_FP in [wgse.install_FP, wgse.prog_FP, "/"] or \
           (wgse.reflib and (dir_FP == wgse.reflib.default_FP or dir_FP == wgse.reflib.FP)) or \
           (wgse.BAM and dir_FP == wgse.BAM.file_FP) or \
           (wgse.outdir and dir_FP == wgse.outdir.FP):
            # We want to make sure they do not specify an important area for the installation.  So we do not wipe
            # out important files when we initialize the temporary directory area. So do not allow the installation
            # directory, reference library, BAM source directory, or output directory areas (if set).
            # Or maybe simply a bad path (not a directory, does not exist, etc)
            # Trying to be extra careful as we delete any files or directories in the Temporary Files directory
            if wgse.lang and wgse.lang.i18n:
                wgse_message("error", 'InvalidTempDirTitle', True,
                             wgse.lang.i18n['errTempDirPath'].replace('{{tempf}}', dir_FP))
            else:
                DEBUG(f'***ERROR: Bad Temporary File Path: {dir_FP}')
            return

        # OK, finally, everything is good and we can set the values
        self.save_FP = dir_FP
        self.set = self.save_FP != self.default_FP
        self.FB  = os.path.basename(dir_FP[:-1]) # Trick; remove trailing slash so basename returns directory name
        # self.list[] is independent of the temporary files directory setting

        if not self.default_FP:     # If just initializsing, default is not set yet. Set to dir without PID yet
            self.default_FP = dir_FP

        # Adjust for process ID (runtime) based temporary directory
        pid = str(wgse.os_pid) + "/"
        if dir_FP + pid == self.FP:     # No change; nothing to do and not initializing. Silently return.
            return
        self.FP  = dir_FP + pid
        self.oFP = nativeOS(self.FP)

        if not os.path.isdir(self.oFP):  # Maybe resetting to old temp area where pid directory already exists
            try:
                os.mkdir(self.oFP)  # Create PID-based temporary directory if not there already
            except:
                if wgse.lang and wgse.lang.i18n:
                    wgse_message("error", 'InvalidTempDirTitle', True,
                             wgse.lang.i18n['errTempDirPath'].replace('{{tempf}}', self.oFP))
                else:
                    DEBUG(f'***ERROR: Cannot create Temporary File directory: {self.oFP}')

        # Want to start with an empty temporary files directory to make sure not the wrong area we wipe out
        # So warn user if not empty.
        dir_contents = [x for x in os.listdir(self.oFP) if not x.startswith('.')]
        if len(dir_contents) > 0:
            if wgse.lang and wgse.lang.i18n:
                wgse_message("warning", 'InvalidTempDirTitle', False, 'errTempDirNotEmpty')
            else:
                DEBUG(f'+++WARNING: New Temporary File area already has content'
                      f'(# files: {len(dir_contents)})\n{dir_contents}')

        self.update_free()

    def update_free(self):
        if self.oFP:
            self.free = shutil.disk_usage(self.oFP).free
        return self.free

    @staticmethod
    def clean_item(oFN, ignore_debug_mode=True):
        """ A simple clean item call with all error check done before.  Only call directly if done so. """

        force = not wgse.DEBUG_MODE and not ignore_debug_mode

        try:
            if os.path.isdir(oFN):      # yleaf tempYleaf dirctory is only commonly occuring directory
                oFN += wgse.os_slash
                DEBUG(f"Clean_Temp dir: {oFN} (DEBUG_MODE = {wgse.DEBUG_MODE})")
                shutil.rmtree(oFN) if force else ''  # Must be careful; hierarchical tree removal
            elif os.path.isfile(oFN):
                DEBUG(f"Clean_Temp file: {oFN} (DEBUG_MODE = {wgse.DEBUG_MODE})")
                os.remove(oFN) if force else ''
            else:
                DEBUG(f"Cleanup element exists? {oFN} (ignoring)")
        except:     # Nobody likes except-all but best to just ignore if any of a number of things happen ..
            pass

    def clean(self, incl_temp_dir=True):        # To clean the list[] and temp dir
        """
        Temp Directory clean-up call.  First the temporary list.  Then the temporary (process) folder itself.
        Do not delete anything if DEBUG_MODE on (just report planned action).
        """

        # First, cleanup files on the clean-up list
        for oFN in self.list:
            # Todo maybe check if file to clean is in install, reflib (new), outdir, or BAM file location trees
            self.clean_item(oFN, ignore_debug_mode=False)
            # wgse_message("error", 'InvalidTempDirTitle', True,
            #              wgse.lang.i18n['errTempElement'].replace('{{ELEM}}', oFN))

        self.list.clear()       # Clear list after cleaning each individual item

        if not incl_temp_dir:
            return

        # Just being extra careful before we do a rm all on a temp directory (process subdirectory!)
        if not (self.oFP and os.path.isdir(self.oFP)):
            wgse_message("error", 'InvalidTempDirTitle', True,
                         wgse.lang.i18n['errTempDirPath'].replace('{{TEMPF}}', self.oFP))
            return

        # Now clean-up the temp directory itself if DEBUG mode is off -- want to save files for analysis if DEBUG mode is on
        for oFBS in os.listdir(self.oFP):
            self.clean_item(self.oFP + oFBS, ignore_debug_mode=False)    # make absolute path before cleaning


######################################################################################################################
class LanguageStrings:
    """
      Internationalization subsystem providing language specific labels, buttons and messages via table lookup.

      Entry Points:
        __init__ (Class creation entry; takes in language.csv file spec which is read in now)
        change_language  (waterfall call from read_language_setting; also from mainWindow language button)
        switch_language  (to set a language or switch mid-stream; forces new mainWindow if a language change)

        Mainly for language translators: there is a reload languages.csv button to recreate this class instance and
         start with change_language_setting again (like at program start). We purposely clear the saved language so the
         language query pop-up happens again as well as to help verify (a new) languages.csv file loaded correctly.

      Note:
        Moved get_lanugage GUI to mainwindow.py button_set_language (keep all GUI stuff there)
        Moved store / read language saved setting to generalized functions in settings module (added more settings)
        This greatly simplified this class.
    """

    def __init__(self, language_oFN):
        """
            Read in file language_oFN as init action. Is a CSV dictionary file of language translations
            wait till actually read_language_setting to then set language and an active wgse.lang.i18n dictionary.
        """

        # Class variables set now (by languages.xlsx file)
        self.language_oFN = language_oFN
        self._lang = {}         # Dictionay containing the languages.xlsx file section of translations
        self.request = {}       # Local language request strings from languages.xlsx
        self.avail_langs = []

        # Class variables to be set later (once language is selected)
        self.language = None    # Keep null; language not set yet; just the universal language translation table
        self.i18n = {}          # only set to particular translation dictionary once a language is selected
        self.selectLanguageWindow = None        # Used to cancel language window popup from here (convenience)

        # Read in the language.csv file to determine languages supported and their translations
        try:
            self._read_languages_xlsx()
        except:
            print(f"***Error reading from system Language Translations file: {self.language_oFN}")
            raise ValueError

    def _read_languages_xlsx(self):
        """
        Read in Languages.xlsx translation spreadsheet and setup translation dictionary for each language.
        Sets self.avail_langs (row 0), self._lang (rest of languages.xlsx) and self.request (language selection strings)
        """
        import openpyxl

        # Read in languages.xlsx and use iterator to process after
        wb = openpyxl.load_workbook(self.language_oFN)
        ws = wb.active
        dim = ws.calculate_dimension()
        DEBUG(f'Language.xlsx Dimension: {dim}')

        rows = ws.iter_rows(values_only=True)  # Iterator for rows of languages.xlsx spreadhseet

        first = next(rows)  # First row is languages available (column headers); special
        num_of_lang = int(first[0])                                 # Cell A1 is number of languages (# columns)
        avail_langs = [first[i + 1] for i in range(num_of_lang)]    # Read language names (row 0) into array
        self._lang = {key: {} for key in avail_langs}               # Create empty translation dictionary per lang

        # Rest of rows are language translations; first column is the key to each languages translation entry
        for row in rows:  # Iterate through rows (language translation entries)
            k = row[0]  # lang translation entry key
            for i in range(0, len(avail_langs)):  # Iterate through languages / columns
                self._lang[avail_langs[i]][k] = row[i + 1].replace("^", "\n") if row[i + 1] else ""

        # Create dictionary of language request strings for get_language pop-up; no reliance on _lang{}{} then
        for i in range(0, len(avail_langs)):
            self.request[avail_langs[i]] = self._lang[avail_langs[i]]["RequestLanguage"]

        self.avail_langs = avail_langs

    def change_language(self, language=None):
        """
        Called after loading settings from user file. If not set or not valid, query user.
        Local import to avoid loop. All GUI pushed to mainwindow module.
        """
        from mainwindow import button_set_language

        if language not in self.avail_langs:      # If requested language not set or not understand
            button_set_language()           # Ask user for language; calls switch_language() directly
        else:
            self.switch_language(language)  # Switch language

    def switch_language(self, language):
        """
        Language Subsytem function to set/switch language; externally called to change language.
        May be called before initial window is setup; called as result of set_language pop-up dialog and buttons
        """
        from mainwindow import mainwindow_reset, mainwindow_resume

        # If got here from the Language Selection button; destroy its pop-up window (ignore errors)
        # We saved the window ID in the class' local variables
        if self.selectLanguageWindow:
            try:
                self.selectLanguageWindow.destroy()
            except:
                pass
            self.selectLanguageWindow = None

        if language == self.language:           # If no change, then simply return
            pass

        elif language not in self.avail_langs:  # If not valid language to switch too ...
            DEBUG(f'*** INTERNAL ERROR: Unknown language: "{language}"; ignoring')

        else:                                   # Everything OK; setup new language translation for selected lang
            self.i18n = self._lang[language]    # Set wgse.lang.i18n to selected language dictionary
            self.language = language            # Save current language in class
            DEBUG(f"New Language: {language}")
            wgse.save_settings()

            if wgse.window and wgse.dnaImage:   # Only reset mainwindow if exists and setup (wgse.dnaImage is set)
                mainwindow_reset()              # Destroys current window and runs mainwindow_setup() again

        if wgse.window and wgse.dnaImage:
            mainwindow_resume()


######################################################################################################################
class OutputDirectory:
    """
    A class to handle the Output Directory setup.  Mostly for changing it in a consistent way from the users button
    as well as the stored setting.  And to move the functionality out of mainwindow so it can focus on UI.
    """

    def __init__(self, new_FP=None):
        self.oFP   = None   # path (user specified)
        self.FP    = None   # path (user specified) (Unix/Universal format)
        self.oFPB  = None   # path with BAM basename (both user specified)
        self.FPB   = None   # path with BAM Basename (both user specified) (Unix/Universal format)
        self.FB    = None   # path (user specified) basename (not BAM but last directory name in path)
        self.user_set = False
        self.free = 0       # free space available
        self.default_outbase = f"WGSEv{wgse.major_version}"

        if new_FP:
            self.change(new_FP)

    def change(self, new_FP, user_set=False):
        """
          Called to change (or initially set) the Output Directory.
        """

        if new_FP is None or len(new_FP) < 3:
            # DEBUG(f"*** WARNING Bad new outdir path: {new_FP}")
            return

        # Assure trailing os_slash; universalOS value so always single forward slash
        new_FP += '/' if new_FP[-1] != '/' else ""
        new_oFP = nativeOS(new_FP)

        if new_FP == self.FP:
            return      # No change; simply return. Especially if stll None

        # If a valid and legal output path; then set the global variables for this new setting
        if is_legal_path(new_FP) and os.path.exists(new_oFP) and os.path.isdir(new_oFP):
            self.FP = new_FP
            self.oFP = new_oFP
            DEBUG(f"Output Path (final): {self.FP}")
            if wgse.BAM:
                self.FPB  = new_FP  + wgse.BAM.file_FB
                self.oFPB = new_oFP + wgse.BAM.file_FB
            # Like the Unix basename() for directories; save last name in path
            self.FB = new_FP.split('/')[-2] if len(new_FP.split('/')) > 1 else new_FP
            self.user_set = user_set
            self.update_free()

    def update_free(self):
        # Note: shutil.disk_usage() will not work on network drives.  They have to be mounted locally.
        # Otherwise, fall back to shell command: df <dir> | tail -1 | awk -F'{print $4}' but use 7 for 4 on MacOS
        if self.oFP:
            self.free = shutil.disk_usage(self.oFP).free
        return self.free

    def new_default(self, path):
        """
          If outdir not yet set by user, then change it to the new default here for a new BAM (or FASTQ).
          Default is a subdirectory "WGSEvN" in the path; where N is the current major tool version. Create if the
          directory does not exist. If you cannot create it, return error so the file can be unloaded.
          Note: If user dir in saved settings, then user_set is true. Otherwise, it would not be be saved.
        """

        if not path or self.user_set:
            return

        upath = universalOS(path)                   # Universal (unix) format
        upath += '/' if upath[-1] != '/' else ""    # assure a trailing slash
        base = os.path.basename(upath[:-1])         # Make sure no trailing '/' before call

        # Idea is if path is already a default outdir 'WGSEv5' then do not create a nested one; just use it
        new_outdir_FP = upath + (f'{self.default_outbase}/' if base != self.default_outbase else "")
        new_outdir_oFP = nativeOS(new_outdir_FP)

        # Only change if not a user_set (forced) dir and the outdir would actually change
        if new_outdir_FP != self.FP:
            try:
                os.makedirs(new_outdir_oFP, exist_ok=True)      # Will make the directory IF it does not exist
            except OSError as e:
                DEBUG(f"*** BAM Default makedir failed: {new_outdir_oFP} with OSError {e.errno}")
                return

            if not (os.path.exists(new_outdir_oFP) and os.path.isdir(new_outdir_oFP)):  # If still not exist, exit
                DEBUG(f"*** BAM Default makedir failed: {new_outdir_oFP}")
                return

            # Replace current outdir with new one
            self.change(new_outdir_FP, user_set=False)

        # Sets part of output File variables that includes BAM base name (if outdir did not change but BAM did)
        # elif self.FP and wgse.BAM:
        #     self.oFPB = self.oFP + wgse.BAM.file_FB
        #     self.FPB = self.FP + wgse.BAM.file_FB

    def clear(self):
        """
          Decided to leave Class always available. So clear it out as if doing an init'; or delete then recreate.
        """
        self.__init__()


######################################################################################################################
class FontTypes:
    """
    Setup in mainwindow_init() after the first call to Tk(), this stores the fonts that will be used for the UI.
    OS platform dependent. Must be dynamically called by each routine using fonts to get the latest setting.
    """

    def __init__(self, toplevel=True):
        from tkinter.font import families      # Historical reason for isolated import

        if wgse.os_plat == "Darwin":
            self.face   = "STIX Two Text" if "STIX Two Text" in families() else "Times New Roman"
            self.facelg = "Arial Black"
            self.basept = 13
        elif wgse.os_plat == "Linux":
            # WSLG with Ubuntu has very limited fonts so use DejaVu Serif and Ubuntu if so; also needs smaller pt sizes
            self.face   = "Times New Roman" if 'Times New Roman' in families() else "DejaVu Serif"
            self.facelg = "Arial Black" if "Arial Black" in families() else "Ubuntu"
            self.basept = 9
        else:   # Windows; nominal
            self.face   = "Times New Roman"
            self.facelg = "Arial Black"
            self.basept = 12

        DEBUG(f'Fonts (init): face {self.face}, size {self.basept}, available {self.face in families()}')

        self.default_face = self.face
        self.default_facelg = self.facelg           # Cannot be changed in UI or stored settings currently
        self.default_basept = self.basept

        self.table = self._newtable()

    def change(self, newface=None, newfacelg=None, newbasept=None):
        from tkinter.font import families, nametofont       # Historical reason for isolated import

        change = False
        if newface is not None and newface != self.face:
            change = True
            self.face = newface if newface in families() else self.default_face
        if newfacelg is not None and newfacelg != self.facelg:
            change = True
            self.facelg = newfacelg if newfacelg in families() else self.default_facelg
        if newbasept is not None and newbasept != self.basept:
            change = True
            self.basept = newbasept if 6 <= newbasept <= 20 else self.default_basept
            nametofont("TkDefaultFont").config(size=self.basept, family=self.face)

        if change:
            DEBUG(f'Fonts (change): face {self.face}, size {self.basept}')
            self.table = self._newtable()

    def _newtable(self):
        from tkinter.font import nametofont       # Historical reason for isolated import

        # Overall default font for system
        default_font = nametofont("TkDefaultFont")
        default_font.config(size=self.basept, family=self.face)
        if wgse.window:
            wgse.window.option_add("*Font", default_font)

        # Font for Combobox
        default_fixed_font = nametofont("TkFixedFont")
        default_fixed_font.config(size=self.basept)
        if wgse.window:
            wgse.window.option_add("*TCombobox*Listbox*Font", default_fixed_font)

        return {  # Ugly hand-coded window GUI.  This at least gets us a little more generality ...
            '12': (self.face, self.basept),      # For version in various result pages (below min size)
            '12b': (self.face, self.basept),     # For version in title marquee of program
            '13': (self.face, self.basept + 1),  # Stats by-chr table; Title manual / exit buttons (below min size)
            '14': (self.face, self.basept + 2),  # Minimum size per Marko (for laptops / Macbook with 1366x768 screen)
            '14b': (self.face, self.basept + 2, "bold"),
            '14o': (self.face, self.basept + 2, "overstrike"),
            '14u': (self.face, self.basept + 2, "underline"),
            '16': (self.face, self.basept + 4),
            '16b': (self.face, self.basept + 4, "bold"),
            '18': (self.face, self.basept + 6),
            '18b': (self.face, self.basept + 6, "bold"),
            '20b': (self.face, self.basept + 8, "bold"),
            '28b': (self.face, self.basept + 16, "bold"),
            '32': (self.facelg, self.basept + 20)  # Main Window Program Title; font face is fixed per platform / init
        }


#####################################################################################################################
# Program Version Check and Manipulate Utilities (class?)
def load_json(file):
    """
    General, common json load. Returns json content as dictionary
    """
    import os
    import json
    import requests

    if file is None or not(file[:4] == "http" or os.path.exists(file)):
        DEBUG(f'*** Cannot find JSON file {file}')
        return {}

    URL = file[:4] == "http"

    try:
        if URL:
            json_dict = json.loads(requests.get(file).text)
        else:
            with open(file, "r") as f:
                json_dict = json.load(f)
    except:  # Exception here if file processing error on reading existing file
        DEBUG(f'*** Error processing JSON file {file}')
        return {}

    return json_dict


def save_json(file, json_dict):
    import json

    try:
        with open(file, 'w') as f:
            json.dump(json_dict, f)
    except:
        DEBUG(f"*** Error writing json file: {file}")


def load_wgse_version():
    """
    Read the release.json and program/program.json file to setup the __version__ string for the program
    """

    release_json = load_json(f'{wgse.install_oFP}/release.json').get('release', {})
    track = release_json.get('track', "Beta")
    major_version = ("v" if track != "Dev" else "") + str(wgse.major_version)  # two v's in a row looks weird

    program_json = load_json(f'{wgse.prog_oFP}/program.json').get('program', {})
    minor_version = program_json.get('version', 0)
    date = program_json.get('date', "1 Jan 1970")

    wgse.manual_url = program_json.get('manual', wgse.manual_url)
    wgse.track = track      # Note: set to stored value now; but may be changed by user before updating the release
    wgse.__version__ = f'{track} {major_version}.{minor_version} ({date})'


def strip_file_URI(path):
    """
    Convert URI for local files into valid path to open. Do nothing if non file protocol (e.g. https://)
    """
    if path is None or len(path) < 9:
        DEBUG(f"Path specification too short: {path}")
        return ""

    elif path[:8] == "file:///":    # implied localhost network name
        if path[9] == ":":          # Windows disk form (e.g. C:) so strip all leading slashes
            return nativeOS(path[8:])
        else:                       # Network name so leave double slash to start
            return nativeOS(path[7:])

    elif path[:7] == "file://":     # Specified a network name
        if path[8] == ':':          # Cannot have a disk specified after a specified network name on Windows
            DEBUG(f"Invalid path specification: {path}")
            return ""
        else:
            return nativeOS(path[5:])

    elif path[:6] == 'file:/':      # Specified a local file path
        if path[7] == ":":          # Windows disk form (e.g. C:) so strip all leading slashes
            return nativeOS(path[6:])
        else:                       # Local Unix style file name so leave leading slash
            return nativeOS(path[5:])

    return path         # Not a file:/ spec so simply return as is


def is_wgse_outdated():
    """
    Will check the indivdual versions of latest against the current to determine if this release is up to date.
    Returns True if up to date, false if older. Checks all 4 (or 6 on Windows) version files.
    """
    # Reflib is special because it could have been relocated
    reflib_oFP = wgse.reflib.oFP if wgse.reflib and wgse.reflib.valid else f'{wgse.install_oFP}/reference/'

    # Find and load latest_json release version file based on release.json track and URLs
    release_json = load_json(f'{wgse.install_oFP}/release.json').get('release', {})
    track = release_json.get('track', 'Beta')
    tempURL = release_json.get('baseURL',
                release_json.get('_commentbaseURL',
                  'https://raw.githubusercontent.com/WGSExtract/WGSExtract-Dev/master/'))
    baseURL = strip_file_URI(tempURL)

    latest_json_file = 'latest-release-' + track + '.json'
    tempURL = release_json.get(track + 'URL', baseURL + latest_json_file)
    latestURL = strip_file_URI(tempURL)

    latest_json = load_json(latestURL)

    # Get all the currently installed release versions
    installer_json = load_json(f'{wgse.install_oFP}/scripts/installer.json').get('installer', {})
    program_json = load_json(f'{wgse.prog_oFP}/program.json').get('program', {})
    tools_json = load_json(f'{wgse.install_oFP}/jartools/tools.json').get('tools', {})
    reflib_json = load_json(f'{reflib_oFP}reflib.json').get('reflib', {})

    cygwin64_json = bioinfo_json = {}
    if wgse.os_plat == "Windows":
        cygwin64_json = load_json(f'{wgse.install_oFP}/cygwin64/cygwin64.json').get('cygwin64', {})
        bioinfo_json = load_json(f'{wgse.install_oFP}/cygwin64/usr/local/bioinfo.json').get('bioinfo', {})

    # Could simply do one big conditional that is in the return statement. But created arrays so we could DEBUG print
    current = [installer_json.get('version', -1),
               program_json.get('version', -1),
               reflib_json.get('version', -1),
               tools_json.get('version', -1)]
    latest = [latest_json.get('installer', {}).get('version', 0),
              latest_json.get('program', {}).get('version', 0),
              latest_json.get('reflib', {}).get('version', 0),
              latest_json.get('tools', {}).get('version', 0)]
    if wgse.os_plat == "Windows":
        current += [cygwin64_json.get('version', -1),
                    bioinfo_json.get('version', -1)]
        latest += [latest_json.get('cygwin64', {}).get('version', 0),
                   latest_json.get('bioinfo', {}).get('version', 0)]

    DEBUG('Release version arrays: installer, program, reflib, tools [, cygwin64, bioinfo]')
    DEBUG(f'Current: {current}')
    DEBUG(f'Latest: {latest}')

    # Check and return status of current version versus latest release available
    return not all(cur >= lat for cur, lat in zip(current, latest))


def update_release_track(new_track):
    """
    Will update release.json file to a new track; allowing the update to proceed with a new track.
    """
    release_json_file = f'{wgse.install_oFP}/release.json'

    full_json = load_json(release_json_file)
    cur_track = full_json.get('release', {}).get("track", None)
    if cur_track is not None and cur_track != new_track:
        full_json['release']["track"] = new_track
        save_json(release_json_file, full_json)


# Moved update_release() to commandprocessor
