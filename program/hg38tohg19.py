#!/usr/bin/env python3
# coding: utf8
#
# hg38tohg19 liftover tool
# Todo work to generalize for all liftover directions and for any vendor file (not just AllSNPs as currently implied)
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#  (stand alone ; but relies on the WGS Extract support utilities and settings)
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.import os

import sys      # for argv[] if called as separate process / program
import os.path  # for os.remove, os.path.exists

from pyliftover import LiftOver     # Main heavyweight lifter

from utilities import DEBUG, nativeOS, universalOS
from commandprocessor import run_bash_script
import settings as wgse     # Imports WGSE environment as if sub-module (to get ref lib liftover chain file)

import string
notabs = string.whitespace.replace("\t", "")

"""###################################################################################################################
  Wrapper for PyLiftover which is itself based on the UCSC Liftover code and liftover reference file.
    Called by module microarray on BAM extracted CombinedKit that was based on reference model Build38.
    Replaces file it is called with (inplace); often an ~50 MB, 2+million SNP file is extracted from a WGS BAM.
    Although currently called as stand-alone, does pull in some WGSE support modules for consistency
    
  UPDATE v5: Generalized to handle microarray source files from vendors. So may see chromosome 23-26, XY, diff sort
"""
# Todo need to inline this code with main microwindow.py file so does not create new environment


# Need to change specialized names from vendors to common names used in liftover tools
def map_chrom(chrom="0"):
    """
        Now that we are generalizing this routine for all RAW microarray files, we need to take care of
        the special versions some companies have: 23,24,25 and 26 for X, Y, XY and MT; respectively. 25 and XY
        are for the PAR region of the X and Y.  Other vendors simply return diploid X in the PAR region. So we will
        assume X coordinates for such specified regions.

        The liftover tool always wants HG (chr) naming for the entries. But all microarray files use GRCh (numeric).
        So beside mapping the odd names given above, we need to make sure the names are mapped to HG.
    """
    map_special = {"23": "X", "24": "Y", "25": "X", "26": "MT", "XY": "X"}  # Assume all XY / 25 are in X coordinates

    special = map_special.get(chrom, "0")     # Check if we need to map the chromosome; if not returns "0"
    if special == "0":                        # unmapped if chromosomes 1-22, X, Y, and MT (mode = "0")
        lchrom, mode = chrom, "0"
    else:                                     # else use the mapped chromosome name; mode is the name before mapping
        lchrom, mode = special, chrom

    lchrom = ("chr" + lchrom).replace("chrMT", "chrM")      # Convert to hg (chr) naming

    if lchrom in wgse.valid_chromosomes:
        return lchrom, mode     # Return liftover loolup chromosome name; and mode (whether changed and from what)
    return "", ""


# Todo generalize to more source / target builds
def my_hg38tohg19(file_oFPB):
    """
    Glue and single call for current use based on HG38 to HG19 liftover of allSNPs TSV microqarray file.
    Overwrites passed in file with new data.  Retains header, chromosome name format, etc.  Just new coordinates.
    """

    # Todo handle .zip file as well as .txt
    #  import zipfile
    #  with ZipFile('aa.zip') as myzip:
    #    with myzip.open('aa.txt') as myfile:   # Use myzip.namelist() to check only 1 file and get correct name?
    # Todo put tmp, hdr and -sorted.txt files into temp/procid folder (available in stand alone mode?)

    # Setup main and intermediate file names (should have been passed the CombinedKit TSV
    txt_oFN        = f'{file_oFPB}.txt'         # Input and eventual output (replaced)
    tmp_oFN        = f'{file_oFPB}.tmp'         # Intermediate; after liftover
    head_oFN       = f'{file_oFPB}.hdr'         # Saved header to prepend after sort
    sorted_txt_oFN = f'{file_oFPB}-sorted.txt'  # Intermediate; after liftover and sort; before rename back to original

    file_FPB       = universalOS(file_oFPB)     # Prepare for quoted, universalOS versions
    txt_qFN        = f'"{file_FPB}.txt"'
    tmp_qFN        = f'"{file_FPB}.tmp"'
    head_qFN       = f'"{file_FPB}.hdr"'
    sorted_txt_qFN = f'"{file_FPB}-sorted.txt"'

    if not os.path.exists(txt_oFN):
        print(f"*** Fatal error. File to liftover does not exist: {txt_oFN}")
        exit()

    # Todo could auto download reference in pyliftover; or by reference library call. In case update available.?
    print("Converting Microarray Build 38 to Build 37/19 (reading the liftover database) ...")

    # Instantiates Liftover class from PyLiftover
    lo = LiftOver(wgse.reflib.liftover_chain_oFN)

    print(f"Finished reading the liftover database; starting the liftover of {txt_qFN}")

    # Do liftover with txt as input, tmp as output
    bad_chrom = 0
    bad_pos = 0
    header = ""

    f_source = open(txt_oFN, "r")       # microarray format file to read in and iterate over

    # Read header lines first ; write out to temp file as need to sort data before tacking back on to head of file
    f_header = open(head_oFN, 'wb')

    source_line = None          # In case the file is empty
    for source_line in f_source:

        line = source_line.strip()
        check_line = line.upper()

        # Silently skip blank lines (shouldn't be any but ...)
        if not check_line:
            pass

        # Only single character line acceptable is blank header line starting with a hash (#)
        elif len(check_line) == 1 and check_line[0] != "#":
            DEBUG(f'Bad source line, skipping: {line}')

        # Check for header line; print to header file and go get next line
        elif check_line and (check_line[0] == "#" or check_line.startswith('"RSID') or
                             check_line.startswith('RSID') or check_line.startswith('SNP')):
            header += source_line
            f_header.write(f'{source_line}'.encode())

        # Start of data area; break out with current source line intact
        else:
            break

    f_header.close()

    # Add a 23andMe v3 header if no header found in the original source file
    if not header:
        head_qFN = f'"{wgse.reflib.cma_FP}raw_file_templates/head/23andMe_V3.txt"'

    split_char = "\t" if "\t" in source_line else "," if '","' in source_line else ""   # TSV or CSV?

    # Read in the rest of file (data lines) and write lifted over coordinate data lines to a temporary file
    f_sink = open(tmp_oFN, "wb")

    while source_line:

        line = source_line.strip(notabs)        # Strip leading / trailing whitespace except tabs
        line_tabs = line.split(split_char)      # split by tabs or commas depending on file type

        # Skip any embedded header, blank line and require atleast four fields (all are ill formatted entries)
        if len(line) < 1 or '#' == line[0] or len(line_tabs) < 4:
            DEBUG(f'Skipping unknown format line: {line}')
            continue

        # Name fields in the split line (if used)
        # origid     = line_tabs[0]
        origchrom  = line_tabs[1]
        origpos    = line_tabs[2]       # Really the only field we will change
        # origresult = line_tabs[3]
        # if len(line_tabs) == 5:         # For two column result files; merge them back together
        #     origresult += f'\t{line_tabs[4]}'

        lchrom, mode = map_chrom(origchrom)     # Map the chromosome name to something liftover understands
        if not lchrom:
            bad_chrom += 1
            DEBUG(f'Skipping invalid chromosome name: {origchrom}, {origpos}')
            continue

        # Actual call to pyliftover to convert positions ; requires UCSC / HGP chrom (chr) name convention
        newcoord = lo.convert_coordinate(lchrom, int(origpos))

        # Successful liftover; write out same line but with new coordinate position
        if newcoord and len(newcoord[0]) > 2:
            # newchrom  = newcoord[0][0]        # todo should check if name is same as lchrom?
            newpos    = newcoord[0][1]
            # newstrand = newcoord[0][2]        # No strand processing now

            new_line = source_line.replace(origpos, str(newpos))
            f_sink.write(f'{new_line}'.encode())

        else:
            bad_pos += 1
            print(f'---WARNING: Cannot liftover position (ignoring): {line}')

        source_line = next(f_source)   # used "for" iterator on file already so cannot use readline()

    f_sink.close()      # Close temporary output data file (lifted over coordinate result)

    f_source.close()    # Close source file

    if bad_chrom or bad_pos:
        print(f'---WARNING: Some entries failed to lift over ({bad_chrom} chromosome names, {bad_pos} positions)')
    else:
        print(f'Succesful completion of liftover for file {txt_qFN}')

    # Going to overwrite the original input file so just delete now before trying to overwrite
    os.remove(txt_oFN) if os.path.exists(txt_oFN) else ''   # Remove initial / final file in case cannot write over

    # Sort the resultant output and prepend with a header before making final file
    # We have to sort as some positions could be translocated
    # Todo some vendor formats are sorted 1, 10, 11, ... instead of 1, 2, 3, 4 ... ;
    #  need to sort within each chromosome to keep orig file chromosome order (especially MT, X, XY, Y)
    command  = (
        f"{wgse.sortx_qFN} -t $'\\t' -k2,3 -V {tmp_qFN} > {sorted_txt_qFN}\n"
        f"{wgse.catx_qFN} {head_qFN} {sorted_txt_qFN} > {txt_qFN}\n"
    )

    run_bash_script("LiftoverCleanup", command)

    # Cannot rely on temp folder for this stand-alone program. So remove directly
    os.remove(tmp_oFN)
    os.remove(sorted_txt_oFN)
    os.remove(head_oFN)

    # Final liftover result is in original source file; simply updated txt_oFN with new coordinate positions


if __name__ == '__main__':
    """ hg38tohg19 liftover MAIN; for when run directly  """

    wgse.init(False)        # Setup local WGSE environment when run stand-alone

    if len(sys.argv) != 2:
        print(f'*** Usage: {sys.argv[0]} <Miroarray RAW File>')
        exit(1)

    # Process command line arguments
    allSNPs_oFPB = nativeOS(sys.argv[1])          # passed as FPB; note in Output and not Temp directory
    # source_ref         = sys.argv[2].strip()

    my_hg38tohg19(allSNPs_oFPB)
