#!/usr/bin/env bash
# WGS Extract v4 in-place common install script (all platforms; was Upgrade_UbuntuLinux.sh in v3)
# Copyright (C) 2021-2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#

# We try to do OS specific stuff in the OS main installer before calling this common installer.
# Serves as a first-time install, upgrade from v2 to v3 to v4, as well as minor release updater (re-entrant)
# Cleans out old release files from v3 and v2 that are no longer used.
# Will fill in any missing reference genomes by user request at the end IF a Fresh install

# We purposely pass in some parameters that we could easily figure out; just to keep user from running by accident
if [[ $# -lt 1 || $# -gt 2 ]] ; then
  printf "Usage: %s arch [maj.min]\n" "$0"
  printf "  arch is either arm64 or x86_64.  maj.min are the major & minor OS release values.\n"
  printf "This should only be called internally from an WGSE OS-specific install script.\n"
  (return 0 2>/dev/null) && return || exit
fi
# echo "Starting $0 $* with $# parameters"
# echo Using Shell `ps -p $$ | grep -v PID | xargs | cut -d" " -f4`

# Common environment setup for scripts here; sets some variables used later so we declare first so shellcheck knows set
declare bashx
declare reflibdir
declare replace
declare upgrade
source scripts/zcommon.sh dummy

# cpuarch=$(uname -m)   # Now passed in as first parameter just so we can require a parameter for this internal script
cpuarch=$1
case $OSTYPE in
  linux*)   cpuarch=$2  ;; # Instead use the Ubuntu major release number for the cpuarch variable
  darwin* | msys* | cygwin*)  ;;
esac

echo
echo '======================================================================================================'
# Setup Python libraries we need (universal except for way to invoke pip)
# On MacOS, we are finding some cases where PIP uses the wrong arch; so have to explicitely set to override
echo "Installing and Upgrading Python libraries on ${OSTYPE}:$1 ... see temp/pip_install.log for details"

opt="--no-warn-script-location"

cmdlist=( install Pillow pyliftover pyscreenshot openpyxl pandas psutil multiqc wakepy )
# bio (brings biopython, numpy, requests urllib3), pyfaidx, pysocks, pySam (python v2 only), elevate

case "${OSTYPE}:${cpuarch}" in    # MacOS passes arch as first arg; maj/min version as 2nd/3rd; Ubuntu maj as first
  darwin*:arm64*)           PIP=( arch -arm64  /usr/local/bin/pip3 )
                              cmdlist+=( tkmacosx "$opt" )  ;;  # MacOS extra package for tkinter colored buttons
  darwin*:x86_64*)          PIP=( arch -x86_64 /usr/local/bin/pip3 )
                              cmdlist+=( tkmacosx "$opt" )  ;;
  linux*:18*)               PIP=( sudo -H python3 -m pip )      # 18.x Ubuntu pip errors out
                              opt=""  ;;                        # $opt not recognized in 18.x Ubuntu pip
  linux*:20* | linux*:22*)  PIP=( pip3 )
                              cmdlist+=(          "$opt" )  ;;
  msys* | cygwin*)          PIP=( python/python.exe -m pip )    # pip3 in cygwin64 requires python/ and
                              cmdlist+=(          "$opt" )  ;;  #  python/scripts be on the path
  *)  printf "*** Error: unknown OS:ARCH combination of %s:%s\n" "$OSTYPE" "$cpuarch"
      (return 0 2>/dev/null) && return 1 || exit 1 ;;
esac

[ ! -e temp ] && mkdir temp     # Make sure temp/ folder is there to take pip_install.log; could put in scripts/ ?

IFS=" " read -r -a cmd <<< "${PIP[@]}" ; cmd+=(install --upgrade pip "$opt")
"${cmd[@]}" | tee    temp/pip_install.log | grep -v "Requirement already satisfied\|^Collecting\|Preparing\|Running\|Using legacy\|Using cached"

IFS=" " read -r -a cmd <<< "${PIP[@]}" ; cmd+=( "${cmdlist[@]}" )
"${cmd[@]}" | tee -a temp/pip_install.log | grep -v "Requirement already satisfied\|^Collecting\|Preparing\|Running\|Using legacy\|Using cached"
# We use the GREP to remove the common, success strings on the many packages and dependencies.  Slims down output.

# Todo MultiQC has bug with Windows release; cannot run from different disk than files it operates on
#   so need to patch the python code in the library to get around / enable.
#   See https://gitter.im/ewels/MultiQC?at=622591f9ddcba117a20d53d4

echo '... finished upgrading Python 3.x libraries'

# NOTE: in v4, pulled the Reference Library out to a seperately installed, versioned ZIP. No longer use a full ZIP.
#  Microarray templates moved into the Reference Library as well. Also moved yLeaf, Haplogrep and FastQC to a
#  seperate tools ZIP to further conserve downloads.  The reference library and tools tend to be more static. Each
#  ZIP has its own version json file to control when it gets updated. Finally, the installer itself is given a
#  version file and special action to restart the installer if it is updated.

echo
echo '======================================================================================================'
# Check and Upgrade WGS Extract installer scripts; if needed
install_or_upgrade installer

echo
echo '======================================================================================================'
# Check and Install or Upgrade main WGS Extract program scripts
run_library=false
install_or_upgrade program
$replace && ! $upgrade && run_library=true

echo
echo '======================================================================================================'
# Check and Install or Upgrade WGS Extract Reference Library templates and scripts files
install_or_upgrade reflib
proc_refgenomes="${WGSEFIN}/scripts/process_refgenomes.sh"  # Moved from reference/genomes/process_reference_genomes.sh

echo
echo '======================================================================================================'
# Check and Install or Upgrade WGS Extract local tools (yLeaf, FastQC, jartools/haplogrep, etc)
install_or_upgrade tools

\rm -f latest.json || true    # Now retained between calls to install_or_upgrade as a common file between packages

echo
echo '======================================================================================================'
echo 'Cleaning up from any previous releases'
# Perform v2 to v3/v4 transformation of any Blobs we want to transfer, recreate or download
# Remove any remnants of the old v2 and v3 install we no longer use

# Handle removing old start / install scripts: 2b release and patches

# Avoiding compiled Applescript due to Translocation issues and Apple not allowing distribution outside signed app
\rm -rf Install_MacOSX.app Start_MacOSX.app Uninstall_MacOSX.app || true
\rm -f Install_MacOSX.scpt Start_MacOSX.scpt Uninstall_MacOSX.scpt || true

# renamed to just MacOS due to BigSur 11 and renamed .sh to .command for easier click-start
\rm -f Install_MacOSX.sh Start_MacOSX.sh Uninstall_MacOSX.sh || true

# Changed all OS specific START files to WGSExtract.xxx and changed from .sh to .command on MacOS
\rm -f Windows_START.bat MacOS_START.sh Linux_START.sh || true

\rm -f 00README.txt "WGSE Betav3 Release Notes.txt" set_WGSEpath.bat  || true
\rm -f Upgrade_v2tov3.command Upgrade_v2tov3.sh Upgrade_v2tov3.bat || true

\rm -f WGSExtractv2b_Francais_Patch.zip WGSExtractv2b_MacOSX_Patchv3.zip WGSExtractv2b_MacOSX_Patchv4.zip || true

# Finished with saving files from the old 2b programs/ folder; we can now remove it
[ -d programs ] && \rm -rf programs # Removes old Win10 binaries as well

# Some old releases delivered corrupted ACLs for temp; so lets try to correct for that
[ -d temp ] && \rm -rf temp
mkdir temp
chmod a+rwx temp

#
# In v2 to v3, change from reference_genomes folder to reference/genomes as default. But we also use the relocated
# reference library location if it exists (should not if v2 to v3 but would if v2 to v4 or later)
#

# Move v2 Reference Genomes to new v3 area (the only reference items saved from v2 release)
# Should be guaranteed reference/genomes in installation. But just in case v4 settings file exists to change it ...
newlib="${reflibdir}/genomes/"      # As opposed to simply reference/genomes
if [[ -d reference_genomes && -d "${newlib}" ]]; then
  echo 'Saving existing reference genomes from v2 release ...'
  cd reference_genomes || echo '***Internal ERROR: cd reference_genomes'
  [ -f hg38.fa.gz ] && [ ! -f "${newlib}hg38.fa.gz" ] && echo Moving hg38 && mv -f hg38.fa.gz "${newlib}"
  [ -f hs37d5.fa.gz ] && [ ! -f "${newlib}hs37d5.fa.gz" ] && echo Moving hs37d5 && mv -f hs37d5.fa.gz "${newlib}"
  [ -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ] && [ ! -f "${newlib}hs38.fa.gz" ] && echo Moving hs38 && mv -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz "${newlib}hs38.fa.gz"
  [ -f human_g1k_v37.fasta.gz ] && [ ! -f "${newlib}human_g1k_v37.fasta.gz" ] && echo Moving human_g1k_v37 && mv -f human_g1k_v37.fasta.gz "${newlib}"
  [ -f hg19.fa.gz ] && [ ! -f "${newlib}hg19_wgse.fa.gz" ] && echo Moving hg19_wgse && mv -f hg19.fa.gz "${newlib}hg19_wgse.fa.gz"
  cd ..
  echo '... saved needed reference files. Removing v2 release reference genomes directory.'
  \rm -rf reference_genomes
  cd "${newlib}" || echo '***Internal ERROR: cd reference/genomes'
  echo 'Processing reference genomes for new reference library format'
  ${bashx} "${proc_refgenomes}" hg38.fa.gz hs37d5.fa.gz hs38.fa.gz human_g1k_v37.fasta.gz hg19_wgse.fa.gz
  cd "${WGSEFIN}"
  echo
fi


# =============================================================================================================
# WGS Extract v3 to v4 upgrade special operations

# FastQC is included in the tools package like yLeaf.  Not separately installed here.

# Rename of Reference Genomes between versions v3, v4 Alpha, etc
cd "${reflibdir}genomes" || echo "*** ERROR: cd ${reflibdir}genomes"
if [ -f "${reflibdir}genomes/chm13_v2.0.fna.gz" ]; then         # Location could be changed in v3 or v4
  echo 'Renaming T2T CHM13 reference genome to new standard'
  mv -f chm13_v2.0.fna.gz chm13v2.0.fa.gz
  rm -f chm13_v2.0* || true
  ${bashx} "${proc_refgenomes}" chm13v2.0.fa.gz
fi
if [ -f "${reflibdir}genomes/hg19_wgse.fa.gz" ] ; then
  echo "Deleting hg19_wgse reference genome as in error and outdated"
  rm -f hg19_wgse* || true
fi
cd "${WGSEFIN}" || echo "*** ERROR: cd ${WGSEFIN}"

# Cleanup any v3 to v4 file changes
(
  cd "${reflibdir}"
  \rm -f genomes/*.sh || true           # Moved to scripts/ or deleted
  \rm -f TruSeq_Exome_TargetedRegions_v1.2_GRCh.bed xgen_plus_spikein.GRCh38.GRCh.bed || true  # renamed
)
\rm -f WGSE_Betav3_Release_Notes.txt 00README_WGSEv3.txt || true     # Was left hanging around in v3 release by accident
\rm -f samtools.exe.stackdump || true                    # Accidently left in an early Alpha v4 release
# \rm -f jartools/GenomeAnalysisTK.jar jartools/picard.jar        # Distributed in v3 but never used. Will leave for now.
\rm -f Start_* Upgrade_* Install_UbuntuLinux.sh Install_Win10.bat || true
\rm -rf win10tools program/microarray || true    # Replaced with cygwin64 and reference/microarray; respectively
\rm -f WGSE_Betav4_Release_Notes.txt || true     # Dropped WGSE_ prefix on 27 June 2022 (mid-Alpha 4l)
\rm -f zcommon.sh zinstall_common.sh zinstall_stage2windows.sh zxterm.sh || true   # moved on 27 June 2022 to scripts/
\rm -f program/version.json reference/version.json jartools/version.json || true   # renamed to $package.json
\rm -f cygwin64/version.json cygwin64/usr/local/version.json || true               # renamed to $package.json
\rm -f make_release.txt || true                                                     # Accidently left in v3 after install
(
  cd scripts
  \rm -f zprocess_refgenomes.sh zget_and_process_refgenomes.sh zcompare_refgenomes.sh || true # removed added z from name
  \rm -f get_and_process_refgenomes.sh || true      # Major change to script; became singular (added zlibrary_common.sh)
  \rm -f make_release.sh make_release.txt || true   # Only needed by developers; not end users
)

# We now remove files for OSs that are not installed here; leaving only one OS set of scripts
case $OSTYPE in
  darwin*)
    \rm -f Install_ubuntu.sh Library.sh WGSExtract.sh scripts/zxterm_ubuntu.sh Uninstall_ubuntu.sh || true
    \rm -f Install_windows.bat Library.bat WGSExtract.bat scripts/zinstall_stage2windows.sh Uninstall_windows.bat || true
    \rm -f Install_linux.sh Library_Linux.sh Terminal_Linux.sh WGSExtract_Linux.sh || true ;;
  linux*)
    \rm -f Install_macos.command Library.command WGSExtract.command Uninstall_macos.command || true
    \rm -f Install_windows.bat Library.bat WGSExtract.bat scripts/zinstall_stage2windows.sh Uninstall_windows.bat || true  ;;
  msys* | cygwin*)
    \rm -f Install_ubuntu.sh Library.sh WGSExtract.sh scripts/zxterm_ubuntu.sh Uninstall_ubuntu.sh || true
    \rm -f Install_macos.command Library.command WGSExtract.command Uninstall_macos.command || true
    \rm -f Install_linux.sh Library_Linux.sh Terminal_Linux.sh WGSExtract_Linux.sh || true ;;
esac

# Call Library* to allow user to add reference genomes IF full, new install; not on upgrade or no change
run_library=false     # As of 4.40, since have get_and_process built into missing_refgenome check; disabling this
if "$run_library"; then
  echo
  echo 'Running the Library manager on a new or upgraded Reference Library.'
  library_mngr="${WGSEFIN}/scripts/zlibrary_common.sh"
  if [[ -e "$library_mngr" ]]; then
    ${bashx} "$library_mngr" dummy
  else
    echo '*** ERROR cannot find the Reference Library scripts/zlibrary_common.sh file.'
  fi
fi

# Return to the OS platform specific master installer that called this scrupt for any final word