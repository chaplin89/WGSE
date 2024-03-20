#!/usr/bin/env python3
# coding: utf8
#
# A(utosomal)conv(ert) module
#   Take the AllSNPs microaray file and a target template to create the target microarray file
# Todo rename to generate_target_microarray
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#  (stand alone)
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2024 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys
import os.path
import platform

os_plat = platform.system()
os_slash = "\\" if os_plat == "Windows" else "/"


# From utilities.py; aconv.py is standalone so not reading in other modules
def nativeOS(path_FP):
    """
    From Universal to native OS specific.  Replace os_slash forward slash to back slash if Win10.
    Also handles drive designations buried in pseudo-Unix paths to proper Win10 system paths.
    Nothing to do on non-Win10 systems.
    """
    if os_plat == "Windows" and path_FP:
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

    return path_oFP


""""###################################################################################################################
    Module aconv (autosomal microarray file converter / creator)
    (stand alone)

    Currently setup as a standalone called module. Arg checking is minimal; called internally only
    Assumes microarray AllSNPs file exists; that consists of all possible SNPs from a BAM / CRAM file (or VCF)
    This module is called for each individual file format to create, Will be extracted from the AllSNPs TSV file
    
    See usage information at the end for calling convetion.
    
    IMPORTANT: Kit names used here must match the targets[] in microwindow.py (wgse.settings not used here)

    The microarray template library may have more than one file for a given target.  This is needed when the order
    of chromosome entries in the target file is different than the order in the AllSNPs file.
"""
# Todo rewrite to drop multi-template files. Do sort by chromosome and then assemble by targets' order. Either generate
#  separate files for each chromosome and assemble in proper order; or store in memory and write once in proper order.
# Todo need to inline this code within the microfiles.py or possibly microwindow.py module . So can make use of targets
#  definition there to describe characteristics (chromosome order, file extension and type, genotype style, etc.)

# Globals used in this module; leftover from original implementation Todo should remove need for globals here
target_prefix = None
target_file_without_suffix = None
target_suffix = None
target = None
target_type_name = None

f_AllSNPs = None     # type:
AllSNPs = None
f_targetfile = None
num_target_files = None
last_called_line_read = None

called_chrom = None
called_pos = None
called_result = None

templ_chrom = None
templ_pos = None
templ_id = None

templates_oFP = None


def get_template_elements(template_line):
    global target
    global templ_chrom, templ_pos, templ_id
    template_line = template_line.replace('"', '')      # Remove double quotes

    # Names must match the targets[]
    if target in ["FTDNA_V1_Affy", "FTDNA_V2",  "FTDNA_V3",  "MyHeritage_V1", "MyHeritage_V2"]:
        line_elements = template_line.split(",")
        templ_chrom = line_elements[1]
        templ_pos = int(line_elements[2])
        templ_id = line_elements[0]
    elif target == "23andMe_SNPs_API":
        line_elements = template_line.split("\t")
        templ_chrom = line_elements[0]
        templ_chrom = templ_chrom.replace("chr", "")
        templ_pos = int(line_elements[1])
        templ_id = line_elements[2].strip()   
    elif target in ["LivDNA_V1",  "LivDNA_V2"] or "CombinedKit" in target or \
            "Ancestry" in target or "23andMe" in target:
        line_elements = template_line.split("\t")
        templ_chrom = line_elements[1]
        templ_pos = int(line_elements[2])
        templ_id = line_elements[0]


def chrconv(chrom_to_convert):
    global target
    if target in ["Ancestry_V1", "Ancestry_V2"]:
        if chrom_to_convert not in ["MT", "M", "X", "XY", "Y"]:
            chrom_to_convert = int(chrom_to_convert)
        if chrom_to_convert == 23 or chrom_to_convert == 25:
            chrom_to_convert = "X"
        elif chrom_to_convert == 24:
            chrom_to_convert = "Y"
        elif chrom_to_convert == 26:
            chrom_to_convert = "MT"

    if chrom_to_convert in ["MT", "M"]:
        return 23
    elif chrom_to_convert in ["X", "XY"]:
        return 24
    elif chrom_to_convert in ["Y"]:
        return 25
    else:
        return int(chrom_to_convert)


def pos1_smaller_than_pos2(chr1, pos1, chr2, pos2):
    cchr1 = chrconv(chr1)  ;  cchr2 = chrconv(chr2)
    if cchr1 < cchr2:
        return True
    elif cchr1 > cchr2:
        return False
    elif cchr1 == cchr2:
        return pos1 < pos2


def pos1_equal_to_pos2(chr1, pos1, chr2, pos2):
    if chrconv(chr1) != chrconv(chr2):
        return False
    else:
        return pos1 == pos2


def pos1_bigger_than_pos2(chr1, pos1, chr2, pos2):
    if chrconv(chr1) > chrconv(chr2):
        return True
    elif chrconv(chr1) < chrconv(chr2):
        return False
    elif chrconv(chr1) == chrconv(chr2):
        return pos1 > pos2


def get_next_called_vars():
    global f_AllSNPs, called_chrom, called_pos, called_result, last_called_line_read

    called_line = ""
    for called_line in f_AllSNPs:
        if '#' in called_line:
            continue                # Skip comment lines
        # print (called_line)
        break
    if called_line:
        line_columns = called_line.split("\t")
        called_chrom = line_columns[1]
        called_pos = int(line_columns[2])
        called_result = line_columns[3].strip()
    else:
        last_called_line_read = True


def write_line_output_file(output_result):
    global f_targetfile, target
    global templ_chrom, templ_pos, templ_id

    if target in ["FTDNA_V1_Affy", "MyHeritage_V2"]:
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif target == "FTDNA_V2":
        if templ_id == "rs5939319":
            f_targetfile.write("RSID,CHROMOSOME,POSITION,RESULT\n".encode())
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif target == "FTDNA_V3":
        output_line = f'{templ_id},{templ_chrom!s},{templ_pos!s},{output_result}'               # comma-sep
    elif target == "23andMe_SNPs_API":
        templ_chrom = templ_chrom.replace("M", "MT")
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result}'            # tab-sep
    elif target in ["LivDNA_V1", "LivDNA_V2"] or "23andMe" in target:
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result}'            # tab-sep
    elif target == "MyHeritage_V1":
        if output_result == "CT":
            output_result = "TC"
        elif output_result == "GT":
            output_result = "TG"          
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif "Ancestry" in target:
        if output_result == "--":
            output_result = "00"
        if output_result == "CT":
            output_result = "TC"
        if output_result == "GT":
            output_result = "TG"
        # if output_result == '.':
        #    output_result = "00"
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result[0]}\t{output_result[1]}'  # tab-sep
    else:
        output_line = ''
    f_targetfile.write(f'{output_line}\n'.encode())


def get_target_suffix(target):
    """ Get the file suffix type for the GLOBAL parameter target. """
    # Note: TellMeGen extension is CSV but internally a TSV

    if any(x in target for x in ["FTDNA", "MyHeritage", "TellMeGen", "meuDNA", "Genera"]):
        return ".csv"
    elif any(x in target for x in ["CombinedKit", "23andMe", "Ancestry", "LivDNA", "MTHFR", "SelfDecode", "1240K", "HO"]):
        return ".txt"
    else:
        return ".err"


def concat_files(target_prefix, target, target_suffix, num_target_files):
    """ Create target file by concatinating the header and one or more tail files. Using global paramters. """

    # Read header
    with open(f'{templates_oFP}head/{target}{target_suffix}', "r") as f_template_head:
        head = f_template_head.read()

    # Write output; starting with header
    targetfile_oFN = f'{target_prefix}_{target}{target_suffix}'
    with open(f'{target_prefix}_{target}{target_suffix}', "wb") as f_targetfile:
        with open(f'{templates_oFP}head/{target}{target_suffix}', "r") as f_template_head:
            f_targetfile.write(f_template_head.read().encode())

        # For each target file, open to read, concat into targetfile
        for i in range(num_target_files):
            i_str = str(i+1)
            partfilename_oFN = f'{target_prefix}_{target}_{i_str}{target_suffix}'
            with open(partfilename_oFN, "r") as f_part:
                f_targetfile.write(f_part.read().encode())
            os.remove(partfilename_oFN)


# Heavily relies on global values passed in and out instead of parameters
def convert_adna(target_prefix, target):
    global f_AllSNPs, f_targetfile
    global called_chrom, called_pos, called_result, last_called_line_read
    global templ_chrom, templ_pos, templ_id
    global AllSNPs, templates_oFP, num_target_files
    global target_suffix

    # Parameters in are all globals; reset before each call
    # 4 simultaneous open files; so not using with-as to avoid hideous indentation
    f_AllSNPs = open(AllSNPs, "r")
    f_targetfile = open(f'{target_prefix}_{target}{target_suffix}', "wb")
    f_template_body = open(f'{templates_oFP}body{os_slash}{target}{target_suffix}', "r")

    if num_target_files == 1:
        with open(f'{templates_oFP}head{os_slash}{target}{target_suffix}', "r") as f_template_head:
            f_targetfile.write(f_template_head.read().encode())

    last_called_line_read = False
    get_next_called_vars()

    for template_line in f_template_body:
        get_template_elements(template_line)    # Side effect: sets global templ_chrom, templ_pos, templ_id
        while pos1_smaller_than_pos2(called_chrom, called_pos, templ_chrom, templ_pos) and not last_called_line_read:
            get_next_called_vars()
        if last_called_line_read or templ_chrom == 0 or templ_pos == 0 or \
           pos1_bigger_than_pos2(called_chrom, called_pos, templ_chrom, templ_pos):
            write_line_output_file("00" if "Ancestry" in target else "--")
            continue
        elif pos1_equal_to_pos2(called_chrom, called_pos, templ_chrom, templ_pos):
            write_line_output_file(called_result)
            continue

    f_template_body.close()
    f_targetfile.close()
    f_AllSNPs.close()


#####################################################################################################################
#   Start of Main
#     When called stand-alone as separate program; which is the use currently

if __name__ == '__main__':


    if len(sys.argv) != 5:
        print(f'usage: aconv.py target source target microarray_refdir')
        print(f'       target is the keyword identifying target target and version to generate.')
        print(f'       source is the full name with path of the AllSNPs file (.txt form) to subset from')
        print(f'       target is the final ouptut path and BAM base name of the intended generated target file')
        print(f'       cma_FP is the full path of the reference/microarray/ folder with the templates')
        exit()

    # Simple command line argument capture; only an internal script. File paths are quoted in shell and still here
    target = sys.argv[1]
    AllSNPs = sys.argv[2].strip('"')
    target_prefix = sys.argv[3].strip('"')
    cma_FP = sys.argv[4].strip('"')
    cma_FP += '/' if cma_FP[-1] != '/' else ''

    # Designed as stand-alone program and called as such.  No tie in to wgse directly. So must recreate some values.

    templates_FP = f'{cma_FP}raw_file_templates/'
    templates_oFP = nativeOS(templates_FP)

    target_suffix = get_target_suffix(target)

    if target in ["23andMe_V4", "23andMe_V5"]:
        num_target_files = 2
    elif target == "Ancestry_V1":
        num_target_files = 4
    elif target == "Ancestry_V2":
        num_target_files = 5
    elif target == "FTDNA_V3":
        num_target_files = 3
    else:
        num_target_files = 1

    if num_target_files == 1:
        convert_adna(target_prefix, target)
    else:
        for i in range(num_target_files):
            i_str = str(i + 1)
            convert_adna(f'{target_prefix}_{i_str}', f'{target}_{i_str}')
        concat_files(target_prefix, target, target_suffix, num_target_files)
