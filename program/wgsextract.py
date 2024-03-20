#!/usr/bin/env python3
# coding: utf8
#
# WGS Extract main (top-level) moduleto process args, activate settings and other modules, and call GUI
# Todo allow non-GUI command line operation (once button functions are removed from mainwindow GUI module)
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#  (stand alone)
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
  Main WGS Extract program module. As can see, is pretty devoid of effort.  Simply gets global system setup via
  settings module. Then calls the mainWindow processing subsystem setup and loop to wait on a user button input.
  Included by all other modules that need settings. Hence kept simple as well with local imports inside wgse_main()
"""
import os
import sys
from argparse import ArgumentParser

sys.path.append(os.path.dirname(os.path.abspath(__file__)))  # Needed for Windows Embedded Python


def get_arguments():
    """
    For stand-alone, command line invocation (or even double click in GUI after associating BAM/CRAM file
    types with program)
    Todo implement this feature; template / fragment now
    Todo processing if standalone file name without other args
    Todo help args need internationalization
    """
    from settings import __version__

    parser = ArgumentParser(prog='WGSExtract',
            epilog=("If no options are specified, an interactive GUI mode is assumed. \n"
                    "BATCH MODE FROM THE COMMAND LINE IS NOT YET IMPLEMENTED."))

    # -h, --help is automatic with ArgumentParser
    parser.add_argument("-v", "--version", action='version', version=f'%(prog)s {__version__}')

    # Mutually exclusive arguments (only one of the three can be used; target file to work on)
    # Not required because can be in stored settings from previous run or simply starting GUI without anything
    mutex = parser.add_mutually_exclusive_group()

    mutex.add_argument("-B", "--bam",
            dest="Bamfile", required=False,
            help="BAM file to load")

    mutex.add_argument("-C", "--cram",
            dest="Cramfile", required=False,
            help="CRAM file to load")

    mutex.add_argument("-F", "--fastq", action='append',
            dest="fastq", required=False, nargs='+',
            help="FASTQ file(s) to load (two if paired-end)", metavar="FASTQ_FILES")

    mutex.add_argument("-V", "--vcf", action='append',
            dest="vcf", required=False, nargs='+',
            help="VCF file(s) to load", metavar="VCF_FILES")

    """    mutex.add_argument("filename",
            dest="file", required=False,
            help="File(s) to load")    # Support OS File Explorer file click with a defined extension """

    # Really a required button for batch mode; but optional because without any args starts interactive mode
    parser.add_argument("-p", "--process_button",
            dest="button", required=False,
            help="Button you wish to process (required for batch mode)", metavar="BUTTON_TO_CLICK")

    align = parser.add_argument_group('required_alignment_button_args')

    align.add_argument("-n", "--new_outfile",
            dest="newfile", required=False,
            help="New BAM/CRAM file to create (with a .bam or .cram extension)",
            metavar="NEW_BAM_CRAM_FILE")

    align.add_argument("-r", "--reference_genome",
            dest="refgenome", required=False,
            help="Reference genome file to align to", metavar="REFERENCE_GENOME_FILE")

    # Settings
    parser.add_argument("-o", "--output_directory",
            dest="outdir", required=False,
            help="Folder for all Output files (override)", metavar="OUTPUT_DIRECTORY")

    parser.add_argument("-l", "--library",
            dest="reflib", required=False,
            help="Reference library directory (override)", metavar="REFERENCE_LIBRARY_DIRECTORY")

    parser.add_argument("-t", "--threads",
            dest="threads", required=False, type=int,
            help="Processor threads to use (override)")

    parser.add_argument("-d", "--debug", action='store_true',
            dest="DebugMode", required=False,
            help="Turn on Debug Mode")

    return parser.parse_args()


def wgse_main(mode):
    """ WGS Extract main program start as independent task. Gives a direct call. """
    import settings as wgse
    from mainwindow import mainwindow_setup

    # print("Starting WGS Extract ...")  # Put in each startup shell script so appears before any errors

    # Start WGSE subsystem with all the main program settings (settings.py)
    #   Used to be a class; now just module to remove an additional layer of naming
    mode_gui = True if mode == "GUI" else False
    wgse.init(interactive=mode_gui)     # Includes all subsystem init calls including mainwindow_init

    if mode_gui:
        # Print explanation to command script window, setup main window, let's go
        print(wgse.lang.i18n["ExplainWhyTerminalWinIsOpen"])
        wgse.window = mainwindow_setup()
        wgse.window.mainloop()  # Go wait for a button click on the main window ...


# ***************MAIN PROGRAM*******************
if __name__ == '__main__':

    # args = get_arguments()      # Place holder; not implemented yet. Intent is to allow BAM / CRAM file in command line.

    # Todo If no args, startup interactive; otherwise batch mode and process args

    wgse_main("GUI")


"""
  To use the global settings in other modules, we need the following code:
    import settings as wgse
  
  Instead of a class Settings in settings.py, we instead use a simple module settings.py with all the global 
   variables defined and initialized at the top level there. Then simply define an init() function that is called
   from the wgsextract main area. If using a class, it requires an extra level of naming as there are no globals.

  Both techniques still require the additional line:
    font = wgse.font
  which exists to get a global variable without any context / module path.  Simply shortens the common name in the GUI.
  Works only if it is a static value. As it is not now, we use local font = wgse.fonts.table instead.
"""
