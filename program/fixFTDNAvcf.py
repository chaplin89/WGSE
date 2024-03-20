#!/usr/bin/env python3
# coding: utf8
#
# Fixes FTDNA BigY VCF files that are not VCF Format 4.x compliant
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#  (stand alone; pipeline stream operation that uses stdin / stdout)
#
# Copyright (C) 2022-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""
Standalone script to fix FTDNA BigY VCF files. FTDNA VCF files include key=value pairs in the Filter column
which is not allowed by the VCF spec; especially not the 4.1 they claim to follow in the header.  So we substitute all
Filter entries with these illegal values with a fixed (allowed) indicator that we also add to the header.  Most of the
values can be found elwewhere in the VCF row entry (DP, Qual, etc) so not really losing anything. But gain being
able to use bcftools on the files (to annotate with yBrowse listing, for example).

See https://gist.github.com/RandyHarr/2cbb8901306891e810d698ab27d212ce for shell script use example to annotate
FTDNA VCF files with yBrowse SNP names and haplogroups.
"""

import sys
import os

# We assume called in a pipeline with bcftools view used to pass in complete VCF on stdin; filtered output on stdout
header = True   ;   format_found = False   ;   info_found = False
for line in sys.stdin:
    if header and line[0] == '#':           # Add FILTER entries to header of newly inserted values
        if not format_found and "##FORMAT=<ID" in line:     # Add Filter definitions before first FORMAT
            os.write(1, b'##FILTER=<ID=PASS,Description="Passed all filters">\n')    # pre-defined; not required
            os.write(1, b'##FILTER=<ID=DP1,Description="Single Read">\n')            # Majority of entries with DP
            os.write(1, b'##FILTER=<ID=DPM,Description="More than 750 reads">\n')    # For heavy PCR dup'ped reads
            os.write(1, b'##FILTER=<ID=Q40,Description="Quality below 40">\n')       # Low Waulity filter
            os.write(1, b'##FILTER=<ID=GTL,Description="?">\n')                      # Unknown paramter; float 0-1?
            format_found = True             # Added new FILTER definitions to the header
        if not info_found and "##INFO=<ID" in line:         # Will add new INFO fields with annotation
            os.write(1, b'##INFO=<ID=HG,Number=1,Type=String,Description="Haplogroup if sample is derived">\n')
            os.write(1, b'##INFO=<ID=ISOGG,Number=1,Type=String,Description="ISOGG Haplogroup if sample is derived">\n')
            info_found = True
        if "#CHROM" in line:                # Last line of header; exit header mode (optimization)
            header = False
        os.write(1, str.encode(line))  # print header line
        continue

    # Simply PASS along any VCF lines with the PASS Filter set
    fields = line.split("\t")           # Split line into column fields
    if "PASS" == fields[6]:
        os.write(1, str.encode(line))  # Most entries are PASS; no change needed
        continue

    # Modify non-PASS filter VCF lines to have new, fixed entries added to header
    entries = fields[6].split(";")      # Split FILTER column entries
    for i in range(len(entries)):            # Process each entry; replacing with fixed entry added to header
        entries[i] = "DP1" if "DP=1"  == entries[i] else \
                     "DPM" if "DP="   in entries[i] else \
                     "GTL" if "GTL="  in entries[i] else \
                     "Q40" if "QUAL=" in entries[i] else ""
    fields[6] = ";".join(entries)       # Recreate FILTER column from new entries
    line = "\t".join(fields)              # Recreate Line with new FILTER column
    os.write(1, str.encode(line))
