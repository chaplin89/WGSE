#!/usr/bin/env bash
#
# Subset of function previously in get_and_process_refgenomes.sh. Simplified to be just download and process one
#  reference genome model as specified on the command line.  Note: must match filename listed in
#  reference/genomes/genomes.csv  Interactive, menu query to support Library command moved to zlibrary_common.sh.
#  Main functionality moved to get_and_process_refgenome() function in zcommon.sh
#
# Todo Change to a Python function and make a Reference Library manager tab in the main program.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
# Copyright (c) 2021-22 Randy Harr
#

if (( $# < 1 || $# > 2 )) ; then
  printf "Usage: %s RefGenomeFile [NIH|EBI]\n" "$0"
  printf "  Downloads and processes the requested reference genome file listed in genomes.csv\n"
  printf "  Allows a preferred source repository to be specified."
  printf "  Normally only called internally by the WGS Extract program.\n"
  (return 0 2>/dev/null) && return || exit
fi

WGSEDIR=$(/usr/bin/dirname "${BASH_SOURCE[0]}")   # Get the script location to determine install directory
WGSEABS=$(cd "$WGSEDIR"; pwd)                     # By cd'ing to it, resolve any aliases and symlinks
if [[ $(/usr/bin/basename "${WGSEABS}") == "scripts" ]]; then
  WGSEFIN=$(/usr/bin/dirname "${WGSEABS}")        # In case in scripts/ subdirectory then move up a level
else
  WGSEFIN="${WGSEABS}"    # Removed escaping embeeded spaces ${WGSEABS/ /\\ }/
fi

declare reflibdir
source "${WGSEFIN}/scripts/zcommon.sh" dummy

find_reflibdir  # Sets reflib with value from settings or the default installation location if not set

cd "${reflibdir}"genomes/

# Read in genomes.csv file and transpose rows to columns; sets one array variable per column
declare -a source finalf     # First pytcode[]; last initf[], gurl[], menopt[], sncnt[], snnam[], descr[] not meeded
read_genomes_file

# The second parameter is the preferred server (NIH vs EBI).  Reference Genomes with an alternate will be marked with
#  a source[] field of NIH-Alt or EBI-Alt. The idea is, if we prefer NIH, then do not select an EBI-Alt model. If EBI,
#  then do no select an NIH-Alt model.  Any other value of source is a unique location model.
notsource="none"
if (( $# == 2 )) ; then
  [[ "$2" == "NIH" ]] && notsource="EBI-Alt"
  [[ "$2" == "EBI" ]] && notsource="NIH-Alt"
fi

# Go through all the file names in genomes.csv looking for this file (from the preferred source if applicable)
for ((j=1 ; j<${#finalf[@]} ; j++)); do
  if [[ "$1" == "${finalf[$j]}" && "$notsource" != "${source[$j]}" ]]; then
    # If we find the match and the preferred source, use the index to get other array variables
    get_and_process_refgenome "$j"
    (return 0 2>/dev/null) && return || exit
  fi
done

echo "***ERROR Reference Genome File $1 not found in genomes.csv"


