#!/usr/bin/env python3
# coding: utf8
#
# Long Process Progress Meter module
#  For long operations, attempts to estimate time / effort remaining based on status messages from program
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#  (stand alone)
#
# Copyright (C) 2021-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""
Standalone script to process BWA Alignment stderr output and provide a progress line instead of tens of thousands of
individual progress report lines.  Keeps the final output command window cleaner while providing updated progress.
Takes the total read segments expected as the second parameter after the log file.  Normally you can pipe stderr from
BWA to this module via "2> >(python progress.py - numsegs)". Reads from first file parameter (stdin if -) and writes
to stdout (stderr if '-'). Provides Percent done and Estimated Time to Completion along with frequencing of updates.
Had to make file parameter as pycharm does not allow stdin when standalone debugging. BWA Indexing non progress lines
are passed as is.

Note: On Windows with Python, must set env TZ="" before calling Python. Otherwise, now() and others return UTC

Todo move from console printout (piped) here to updating python progress bar in PleaseWait
"""

import sys
import os
from datetime import datetime

dangling_line = False   # True if last line ended with '\r'; so know to send '\n' before when not a status update
total_segs = 0          # Passed in total segments expected to align (in BWA Align)


def time_stamp():
    raw = datetime.now()
    hours = raw.hour  ;  mins = raw.minute
    # time.strftime and time.strptime does not quite do what we want
    stamp = f"{hours-12}:{mins:02n} PM" if hours > 12 else \
            f"{hours}:{mins:02n} AM"    if 0 < hours < 12 else \
            f"{hours:02n}:{mins:02n} AM"    # Must be hours == 0

    return stamp, str(raw)      # hr:min stamp, full date and time stamp


def time_string(time_left):
    """
    Routine to compose custom days and time (left) string from parameter of total seconds.
    Key is the integer division (//) that truncates the value at each stage.
    """
    fover = "time OVER is" if time_left < 0 else "time left is"

    time_left = abs(time_left)
    days  = time_left // (3600 * 24) ; time_left -= days  * 3600 * 24
    hours = time_left //  3600       ; time_left -= hours * 3600
    mins  = time_left //  60         ; time_left -= mins  * 60
    # secs  = time_left         # No longer included (to fine a granularity)

    days_str = f'        '         if days == 0 else \
               f'{days:>2n} day ,' if days == 1 else \
               f'{days:>2n} days,'
    time_str = f'   {mins:02n} mins'          if hours == 0 else \
               f'{hours:02n}:{mins:02n} hour' if hours == 1 else \
               f'{hours:02n}:{mins:02n} hours'     # dropped {secs:02n}

    return f'{fover} {days_str} {time_str}'


def print_mesg(mesg, progress=True):
    """ When piped input from stdin, write to stderr.  Else if reading in a logfile, then write to stdout
        If a progress message in piped mode, then print with '\r' so overwrite last message.
        If dangling and not a progress message, first print a newline to terminate the last progress message
    """
    global piped, dangling_line

    # Force a linefeed on any outstanding line that was a progress line with '\r' (usually only piped also)
    if dangling_line and not progress:
        mesg = '\n' + mesg
        dangling_line = False

    if piped:
        end = '\r' if progress else '\n'        # progress lines out to console end in '\r' to overwrite the previous
        print(mesg, end=end, file=sys.stderr, flush=True)
        dangling_line = progress

    else:
        print(mesg)                 # all nessages to logfile end in newline ('\n'); even progress messages
        dangling_line = False       # should not be set; just for robustness


def print_align_mesg(total_iter, time_so_far, total_seqs_read):
    """
    Print the current BWA Align progress line based on time and read sequences so far.
    If log_file is dash ('-', stdin), then print update progress line to stderr.  Else dump each line to stdout.

    Each iteration is usually around 1 million read segments. time for an update message is heavily dependent on the
    platform and how fast it can process 1 million read segments.
    """
    global piped, total_segs

    fraction_done = total_seqs_read / total_segs if total_segs > 0 else 0

    total_time = int(time_so_far * total_segs / (total_seqs_read + 0.0001))
    ftime = time_string(total_time - time_so_far)       # time_string() adds text of status (left-to-finish or run-over)

    iter_str = f'{total_iter:>3}' if total_iter < 1000 else f'{total_iter}'

    update = time_so_far / total_iter if total_iter > 0 else 600     # Initial guess of 600 seconds before first update

    stamp = f'({iter_str} iters done, one iter every {int(update):>3} secs, last updateed {time_stamp()[0]})'
    mesg = f'*** {fraction_done:>7.2%} aligned; {ftime} {stamp}'
    print_mesg(mesg, progress=True)
    # Purposely use \r and not (default) newline \n so we overwrite the same progress line when piped
    # But during debug or non-piped (file input instead of stdin); use newline with print and '\n' as not overwriting


if __name__ == '__main__':
    module = os.path.basename(sys.argv[0])
    if len(sys.argv) < 3 or len(sys.argv) > 4:  # Remember first arg is always argv[0] - the module name
        print(f'***ERROR: Exactly two parameters needed for {module} call.', file=sys.stderr, flush=True)
        print(f'   python3 progress.py log_file total_read_segments_to_process', file=sys.stderr, flush=True)
        print(f'     log_file can be "-" to read stdin stream', file=sys.stderr, flush=True)
        print(f'    Optional fourth DEBUG parameter allowed also.', file=sys.stderr, flush=True)
        exit(1)
    log_file = sys.argv[1]          # Log file or '-' if stdin (piped)
    total_segs = int(sys.argv[2])   # passed in total_segs (if known, otherwise is a guess)
    DEBUG = len(sys.argv) == 4      # Any 3rd paramter; could be DEBUG like expected

    piped = log_file == "-"         # if piped from stdin, will print to unbuffered stderr instead of stdout

    dangling_line = False
    print_mesg(f'*** Starting progress meter', progress=False)

    # How to distinguish BWA Index and BWA Align; adjust initial print depending on type
    cur_iter = 0    ;  cur_char = 0     ;  total_char = 0           # BWA Index
    total_iter = 0  ;  time_so_far = 0  ;  total_seqs_read = 0      # BWA Align
    if total_segs > 0:      # Indicates an alignment
        # Print initial status line of zeros as it may take awhile to get the first line
        print_align_mesg(total_iter, time_so_far, total_seqs_read)

    f = sys.stdin if piped else open(log_file, "r")
    for line in f:
        # Receive the stderr as stdin; stderr is unbuffered so flushed; so this looping works per line
        line = line[:-1] if len(line) > 0 and line[-1] == '\n' else line    # Strip trailing newline

        if DEBUG and not line.startswith("[M::"):        # Pass through every line if in DEBUG_MODE
            print_mesg(line, progress=False)

        # Process (sink) lines

        # BWA Align progress line to process; every iteration of aprox 1 million segments
        if line.startswith("[M::mem_process_seqs]"):
            # [M::mem_process_seqs] Processed 1077058 reads in 799.447 CPU sec, 50.357 real sec
            fields = line.split(" ")
            total_iter += 1
            time_so_far += int(round(float(fields[8]), 0))
            total_seqs_read += int(fields[2])
            print_align_mesg(total_iter, time_so_far, total_seqs_read)

        # BWA Index progress line; comes every 10 iterations
        elif line.startswith("[BWTIncConstructFromPacked]"):
            # [BWTIncConstructFromPacked] 10 iterations done. 99999994 characters processed. -- most are 10^7 "chars"
            cur_iter = int(line.split(" ")[1])         # Not used ; only first 2/3 are 10^7 per iter
            cur_char = int(line.split(" ")[4])
            fraction_indexed = cur_char / total_char
            msg = f'*** [BWA Index]{fraction_indexed:>7.2%} indexed (last updated {time_stamp()[0]})'
            print_mesg(msg, progress=True)

        # BWA Index initial data line to understand total
        elif line.startswith("[BWTIncCreate]"):
            # [BWTIncCreate] textLength=6199845082, availableWord=448243540
            total_char = int(line.split(" ")[1].split("=")[1][:-1])  # "textlength" is called "chars" in status lines
            print_mesg(line, progress=False) if not DEBUG else ""   # Pass through

        # Toss these alignment lines (just way too verbose; not information added)
        elif line.startswith("[M::"):
            pass

        # Pass through any line not processed (and if not already printed due to DEBUG mode)
        elif not DEBUG:
            print_mesg(line, progress=False)

    if not piped:
        f.close()

    if dangling_line:
        print_mesg(" ", progress=False)
