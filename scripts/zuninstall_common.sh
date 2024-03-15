#!/usr/bin/env bash
#
# WGS Extract Uninstall Script for All (common)
# Copyright (C) 2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# Note: will be in old BASH as exists in MacOS as MacPorts has been deleted by time this is called

# Common environment setup for scripts here; sets some variables used later
if [[ $# -ne 1 ]] ; then
  printf "Usage: %s reflibdir\n" "$0"
  printf "  reflibdir is the WGS Extract set Reference Library location (possibly redirected by user setting).\n"
  printf "This should only be called from WGSE Uninstaller internal scripts.\n"
  exit
fi

declare home
source scripts/zcommon.sh dummy

reflibdir="$1"    # Had to pass in reflibdir as program jq.exe used in zcommon is gone by now (if MacOS or Ubuntu)

echo '======================================================================================================'
echo 'Now removing the WGS Extract program and reference library itself.'
echo
save_reflib=true
if [ -d "$reflibdir" ]; then
  enclose=$(dirname "$reflibdir")
  refbase=$(basename "$reflibdir")
  read -p "Do you want to remove the WGS Extract Reference Library ($refbase) and its genomes [y/N]? " -n 1 -r ; echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ "$refbase" != "/" && "$enclose" != "$refbase" ]]; then  # Just make extra sure we are not wiping whole disk
      read -p "Really remove the WGS Extract Reference Library ($refbase) [y/N]? " -n 1 -r ; echo
      if [[ $REPLY =~ ^[Yy]$ ]]; then
        cd "$enclose"
        rm -rf "$refbase"
        save_reflib=false
      else
        echo " ... Leaving the Reference Library ${refbase} in place."
      fi
    else
      echo '*** Internal error: Reference Library appears as root directory (/)'
    fi
  else
    echo " ... Leaving the Reference Library ${refbase} in place."
  fi
else
  save_reflib=false
fi

remwgse=false
if [ -d "$WGSEFIN" ]; then
  enclose=$(dirname "$WGSEFIN")
  wgsebase=$(basename "$WGSEFIN")
  read -p "Do you want to remove the WGS Extract program ($wgsebase) and its settings [y/N]? " -n 1 -r ; echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ "$wgsebase" != "/" && "$enclose" != "$wgsebase" ]]; then  # Just make extra sure we are not wiping whole disk
      read -p "Really remove the WGS Extract program (${wgsebase}) [y/N]? " -n 1 -r ; echo
      if [[ $REPLY =~ ^[Yy]$ ]]; then
        cd "$enclose"
        remwgse=true
        if [[ "$save_reflib" == true && "$reflibdir" == "${WGSEFIN}/reference/" ]]; then
          echo 'Saving the Reference Library that is inside the WGS Extract installation by moving it up one level.'
          mv "${wgsebase}/reference" "${wgsebase}_reference"
        fi
        rm -f "${home}/.wgsextract" "${home}/.wgsedebug" "${home}/.wgsewslbwa" || true
      else
        echo " ... Leaving the WGS Extract program folder ${wgsebase} in place/"
      fi
    else
      echo '*** Internal error: WGS Extract installation appears as root directory (/)'
    fi
  else
    echo " ... Leaving the WGS Extract program folder ${wgsebase} in place/"
  fi
fi

echo
echo '======================================================================================================'
echo 'Finished uninstalling WGS Extract and its programs. '
echo
echo 'If you specified a Temp directory outside the WGS Extract install folder,'
echo 'you will need to remove that folder yourself. The WGS Extract program folder itself'
echo 'will not delete (if so requested) until you close this window. Due to race conditions,'
echo 'there may be a few files left that you will have to delete manually.'
echo
# Unlike installer, we will pause here as the script that called us will no longer be there
read -n1 -r -p 'Press any key to close this window (after first scrolling up to review for errors) ...'
{   # Critical code section; force pre-read as this script file may be deleted
if "$remwgse"; then
  # Windows is special because we are trying to delete the cygwin64 installation with bash while
  #  in this bash script. So drop out to a cmd.exe forked process to do the deletion after we leave
  case $OSTYPE in
    darwin* | linux*)
      nohup rm -rf "$wgsebase" &>/dev/null & ;;
    msys* | cygwin*)
      wgsebasedos=$(cygpath -w "$wgsebase")
      # delcmd="start /b rmdir ""${wgsebasedos}"" /s/q"
      nohup cmd.exe /d /c start /b rmdir "${wgsebasedos}" /s/q &>/dev/null  ;;
  esac
fi
# Need to force an immediate exit as rest of the file may no longer exist to read
exit
}