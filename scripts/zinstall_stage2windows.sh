#!/usr/bin/env bash
# WGS Extract v4 Win10 install script (after Cygwin64 installed via Install_windows.bat)
#
# Copyright (C) 2021-2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#

# Because we want to spend minimal effort in CMD.EXE BAT file, we have most of the Windows install
# happening here in this BASH script.  Called once a base CygWin64 environment is set up.

# Todo Originally Cygwin64 Python was too old.  They seem to have caught up now. Consider simply
#  installing python and java with the cygwin release.

# Check for required parameter (mainly to verify called internally and not directly)
if [[ $# -ne 1 ]] ; then
  printf "Usage: %s dummy\n" "$0"
  printf "  dummy is a Windows variable passed in.\n"
  printf "This script should only be called internally from the Windows install script.\n"
  exit
fi

# Common environment setup for scripts here; sets some variables used later so declare for scope (shellcheck)
declare curlx
declare osrelease
declare cpu_arch
source scripts/zcommon.sh dummy

# In case called directly, make sure the path to cygwin64 is setup.  The Windows Install script intalls cygwin64.
case $OSTYPE in
  darwin* | linux*)   echo '*** ERROR: This is an installer for MS Windows systems only!' && exit 1  ;;
  msys* | cygwin*)
    if [ ! -d cygwin64/ ]; then
      echo '*** ERROR: Cannot find the Windows Cygwin64 tools previously installed.' && exit 1
    fi  ;;
esac

echo
echo '======================================================================================================'
# WinPython standalone Python release
# Was simply redistributed to everyone in v2; now only install with Win10 users from source release
if [ ! -d python ]; then
  echo 'Installing Python 3.10.2 ...'

  # https://github.com/winpython/winpython/releases/download/4.6.20220116/Winpython64-3.9.10.0dot.exe  39100
  ${curlx} -o Winpython64.exe "https://github.com/winpython/winpython/releases/download/4.6.20220116/Winpython64-3.10.2.0dot.exe"

  if [ -e Winpython64.exe ]; then
    chmod a+x Winpython64.exe && ./Winpython64.exe -y && rm -f Winpython64.exe    # 7Zip self-extracting archive
    # Only need the Python interpreter; not all the IDE that comes with it
    cd WPy64-31020  &&  mv python-3.10.2.amd64 ../python  &&  cd ..  &&  rm -rf WPy64-31020
    echo "... finished installing Python 3.10.2"
  else
    echo "*** ERROR - failed to download Python 3 release."
  fi
else
  echo 'Python already installed.'
fi

echo
echo '======================================================================================================'
# We check to see if Java command is already available; or if we have locally installed it already
# GATK3 / Picard / VariantQC require JDK8; GATK4, FASTQC, etc require openJDK11+. openJDK17 LTS came out in 2021.
# OpenJDK stopped delivering binaries; switched to adoptium. Maybe should go to Azul like for Mac due to M1.
# Now first number updated with every release (e.g. 17) instead of the second (1.8.xx)

# See if java already on the PATH; if so, determine the version so we can avoid installing it and simply use
unset jre
# Returns path to command if it exists
jre=$(command -v java &>/dev/null)
unset jre8
unset jre17
if [[ -v $jre ]] ; then
  IFS='.'
  # Check the version: take first line, 3rd grouping, strip double quote and dots to read into array of numerals
  read -ra ver <<<$($jre -version 2>&1 >/dev/null | head -1 | cut -d" " -f3 | sed -e 's/"//g')
  if [ "${ver[0]}" -eq 1 ] && [ "${ver[1]}" -eq 8 ]; then
    jre8=$jre
  elif [ "${ver[0]}" -ge 11 ]; then
    jre17=$jre
  fi
fi

# Check if ver 11+ installed; if not install jre17 now
if [ -z "$jre17" ] && ! [ -f jre17/bin/java.exe ]; then
  echo 'Installing Java JRE v17'

  # We do not need a standalone Java release but easier to setup that way on Windows from a script
  openjdk="https://github.com/adoptium/temurin17-binaries/releases/download/"
  ${curlx} -o jre17.zip "${openjdk}/jdk-17.0.2%2B8/OpenJDK17U-jre_x64_windows_hotspot_17.0.2_8.zip"

  if [ -e jre17.zip ]; then
    powershell Expand-Archive -LiteralPath "jre17.zip" -DestinationPath "." -Force && rm -f jre17.zip
    mv jdk-17.0.2+8-jre jre17
    chmod a+x jre17/bin/*.exe jre17/bin/*.dll
    # Program does two similar checks; either simply available or if Win10, in the jre subdirectory
    echo finished installing Java JRE v17
  else
    echo "*** ERROR - failed to download Java JRE 17 release from Adptium."
  fi
else
  echo 'Java v17 already installed.'
fi

# Check if ver 8 installed; if not install jre8 now
if [ -z "$jre8" ] && ! [ -f jre8/bin/java.exe ] ; then
  echo 'Installing Java JRE v8'

  # We do not need a standalone Java release but easier to setup that way on Windows from a script
  openjdk="https://github.com/adoptium/temurin8-binaries/releases/download/"
  ${curlx} -o jre8.zip "${openjdk}/jdk8u345-b01/OpenJDK8U-jre_x64_windows_hotspot_8u345b01.zip"

  if [ -e jre8.zip ]; then
    powershell Expand-Archive -LiteralPath "jre8.zip" -DestinationPath "." -Force && rm -f jre8.zip
    mv jdk8u345-b01-jre jre8
    chmod a+x jre8/bin/*.exe jre8/bin/*.dll
    # Program does two similar checks; either simply available or if Win10, in the jre subdirectory
    echo finished installing Java JRE v8
  else
    echo "*** ERROR - failed to download Java JRE 8 release from Adoptium."
  fi
else
  echo 'Java v8 already installed.'
fi

echo
echo "======================================================================================================"
# Grab our own bioinformatics release to add to the existing cygwin64 Unix tools already installed
# Was simply redestributed in v2 to everyone; in v3 only installed with Win10 users. Was mixed in with /bin
# in v3, now installed in /usr/local so can be easily replaced
install_or_upgrade bioinfo      # leaves the latest.json file for use in later package installs

echo
echo '======================================================================================================'
echo 'Calling Common script to finish WGS Extract install ...'
cd "${WGSEFIN}"
/bin/bash scripts/zinstall_common.sh "${cpu_arch}" "${osrelease}"
# Propogate the exit code back to the Windows Batch file that called this script
exit $?

