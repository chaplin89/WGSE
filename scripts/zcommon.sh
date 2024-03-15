#!/usr/bin/env bash
# WGS Extract v4 in-place common scipt startup (all platforms)
# Copyright (C) 2021-2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#
if [[ $# -ne 1 ]] ; then
  printf "Usage: source %s dummy\n" "$0"
  printf "  dummy is anything; required to avoid direct user click.\n"
  printf "  This should only be called from the WGSE internal scripts.\n"
  (return 0 2>/dev/null) && return || exit
fi

if [ ! "${BASH}" ]; then
  echo '*** ERROR: '"$0"' script must be started from within BASH.'
  (return 0 2>/dev/null) && return || exit
fi

export curlx="curl -kLC - --retry 5"
export zipdir="WGSExtractv4"   # Directory stored in zip files (used in make_release.sh and installer)
export release=WGSExtractv4    # Set toplevel folder inside .zip

# Most of the declare and export statements are not needed as this file is sourced. We add them to document what is
# being set and used externally by scripts that source this one (functions, variables) and to help shellcheck
# There could be instances of sub-shell calls that need them though

# ------------------------------------------------------------------------------------------------
# PATH not fully set (yet); override builtin (ancient) BASH on path for some OSs. Need path for dirname, basename, etc
case $OSTYPE in
  msys* | cygwin*)
    bashx="/bin/bash.exe"
    if [[ ":$PATH:" != *":/usr/local/bin:"* ]]; then
      PATH="/usr/local/bin:/bin:${PATH}"
    fi   ;;
  darwin*)
    bashx="/opt/local/bin/bash"
    if [[ ":$PATH:" != *":/opt/local/bin:"* ]]; then
      PATH="/opt/local/bin:/opt/local/sbin:${PATH}"
    fi   ;;
  linux*)
    bashx="/usr/bin/bash"  ;;
    # All tools installed in system release directories already on path by apt or our local installer
  *)
    printf "*** Error: unknown OSTYPE of %s\n" "$OSTYPE"
    (return 0 2>/dev/null) && return || exit  ;;
esac
export bashx
export PATH         # System wide (BASH) variable

# Helpful advice on distinguishing unset, null and non-null values in a BASH variable
# https://www.geeksforgeeks.org/bash-scripting-how-to-check-if-variable-is-set/
# https://stackoverflow.com/questions/3601515/how-to-check-if-a-variable-is-set-in-bash/16753536#16753536
# https://stackoverflow.com/questions/35006457/choosing-between-0-and-bash-source

# ------------------------SETUP PATH TO INSTALLATION FIRST------------------------------------------
# Make sure we are in the installation directory (and leave the variable there to find it again) ...
if [[ -z ${WGSEFIN:+set} ]]; then           # If not already set ...
  WGSEDIR=$(dirname "${BASH_SOURCE[0]}")    # Get the (calling) script location to determine install directory
  WGSEABS=$(cd "$WGSEDIR"; pwd)     # By cd'ing to it, resolve any aliases and symlinks (readlink not always available)
else
  WGSEABS="${WGSEFIN}"
fi
if [[ $(basename "${WGSEABS}") == "scripts" ]]; then
  WGSEABS=$(dirname "${WGSEABS}")   # In case in scripts/ subdirectory then move up a level
fi
case $OSTYPE in
  msys* | cygwin*)
    temp=$(cygpath -u "${WGSEABS}") # pwd in Cygwin64 sometimes returns DOS style path with trailing slash
    WGSEFIN="${temp%/}" ;;
  darwin* | linux*)
    WGSEFIN="${WGSEABS%/}"  ;;      # REMOVED Escape embedded spaces: ${WGSEABS/ /\\ }
  *)  printf "*** Error: unknown OSTYPE of %s\n" "$OSTYPE"
      (return 0 2>/dev/null) && return || exit  ;;
esac
# cd "${WGSEFIN}"     # cd will accept "sp ace", sp\ ace, but not "sp\ ace"; which we had at first
export WGSEFIN


# ------------------------------------------------------------------------------------------------
# Get users home directory; if in Windows, their true user home and not cygwin BASH one
case $OSTYPE in
  msys* | cygwin*)
    home=$(cygpath.exe -u "$USERPROFILE") ;;     # Below is cmd.exe BAT version
    # %CYGBIN%\cygpath.exe -u %USERPROFILE% > temp.txt  &  set /p home=<temp.txt  &  del temp.txt
  darwin* | linux*)
    home=~  ;;                                  # Autoexpands to the full path when assigned
esac
export home

# ------------------------------------------------------------------------------------------------
# Get OS Release version info (1-3 parts usually; in array)
case $OSTYPE in
  msys* | cygwin*)
    # In Cygwin, the main cygwin DLL is the Unix kernel version.  base-cygwin is more a cygwin release version.
    # Presume we do not want a windows version which is pretty useless anyway; builds are more useful there.
    osrelease=$(cygcheck -c base-cygwin| grep base-cygwin | tr -s "[:space:]" | cut -d" " -f2)
    osrelease=${osrelease//-/.} ;;
  darwin*)
    # Pre BigSur, versions were 10.maj.min.  Since BigSur, version are maj.min with BigSur starting at 11.
    osrelease=$(sw_vers -productVersion)  ;;
  linux*)
    # This only works for Ubuntu that has the lsb-release file
    [ -e /etc/lsb-release ] && IFS='=' read -ra elements <<<$(grep "ION=\"Ubuntu [12]" /etc/lsb-release)
    [ ${#elements[@]} -eq 2 ] && IFS=' ' read -ra parts <<<"${elements[1]}"
    [ ${#parts[@]} -gt 1 ] && osrelease=${parts[1]}
    case ${osrelease} in
      18*|20*|22*)  ;;
      *)  echo '***Warning: We only handle Ubuntu LTS releases: 18, 20, 22; setting to 18'
          osrelease="18.04.0"
    esac  ;;
esac
cpu_arch=$(uname -m)    # Generally, always x86_64. But may be aarm64 for Apple M1/M2's
export osrelease cpu_arch

# ------------------------------------------------------------------------------------------------
# Get Python access (assuming only used when already installed; as this file is sourced by the installer)

case $OSTYPE in
  msys* | cygwin*)
    pythonx="${WGSEFIN}/python/python.exe"
    oWGSEFIN=$(cygpath -m "${WGSEFIN}")  ;;
  darwin*)
    pythonx="/usr/local/bin/python3"
    oWGSEFIN="$WGSEFIN"  ;;
  linux*)
    pythonx="python3"
    oWGSEFIN="$WGSEFIN"  ;;
esac
export pythonx  oWGSEFIN

export process_refgenomes="${oWGSEFIN}/scripts/process_refgenomes.sh"

# ------------------------------------------------------------------------------------------------
#
# Globally used function calls
#
# Any "returned" values are declared as export and defaults set before the function definition
#

# Need to check wgse settings in case the user moved the reference library; relies on other settings in this file
# Needed in uninstall_windows.bat before deleting jq there so replicated this code in cmd.exe style there also
export reflibdir=""
find_reflibdir () {
  reflibdir="${WGSEFIN}/reference/"       # Default location in installation directory

  if [[ -e "${home}/.wgsextract" ]]; then
    # In case called before jq installed (during WGSE installation; by Library command before full install, etc)
    if command -v jq &> /dev/null ; then      # Should have been installed already ...
      newreflib="null"                        # If command below fails; then returns "null", not ""
      newreflib=$(jq -r '."reflib.FP"' "${home}/.wgsextract" )   # Return string from settings (else "null")

      if [[ "$newreflib" != "null" ]]; then     # Because jq returns a 4 character string "null" and not ""
        case $OSTYPE in
          msys* | cygwin*)    tempref=$(cygpath.exe -u "$newreflib")  ;;    # Massage Windows version
          darwin* | linux*)   tempref="$newreflib"  ;;
        esac
        # printf "Reference Library was at %s\n but moved to %s in settings.\n" "${reflibdir}" "${tempref}"
        reflibdir="${tempref}"              # From settings file -- reference library was moved
      fi
    else
      # See https://github.com/stedolan/jq/ for more information.
      echo "*** ERROR: JQ is needed by find_reflibdir but not found."
    fi
  fi

  # echo "find_reflibdir: $reflibdir"
}
export -f find_reflibdir

# ------------------------------------------------------------------------------------------------
# Version check system
#  Similar function needed in install_windows.bat so replicated in cmd.exe there
#  Could potentially be local to zinstall_common.sh but zinstall_stage2windows.sh needs also

declare release_track base_url latest_package_url
get_latest_json() {     # $1 is package; may get and use $base_url
  # We try to be a little smart.  If a latest.json file already there, then check the latest_package_url to see
  # if it changes here.  If so, replace the latest.json with the new one indicated.  If not, simple return.
  saved_package_url="$latest_package_url"   # Save value of current latest.json file; if exists / set

  if [[ -z "$release_track" ]]; then        # Read release.json content for the first time
    release_json="${WGSEFIN}/release.json"

    if command -v jq &> /dev/null && [ -e "$release_json" ] ; then
      release_track=$(jq -r .\"release\".\"track\" "$release_json")
      base_url=$(jq -r .\"release\".\"baseURL\" "$release_json")
      latest_package_url=$(jq -r .\"release\".\""${release_track}"URL\" "$release_json")
    fi
  fi

# Currently just getting Developer track info from combined latest-release file like other tracks
#  if [ "$release_track" == "Dev" && -n $base_URL]; then    # Override for Developer track; from github file system
#    case $1 in
#      installer)  latest_package_url="${base_url}/scripts/installer.json"          ;;
#      program)    latest_package_url="${base_url}/program/program.json"            ;;
#      reflib)     latest_package_url="${base_url}/reference/reflib.json"           ;;
#      tools)      latest_package_url="${base_url}/jartools/tools.json"             ;;
#      bioinfo)    latest_package_url="${base_url}/cygwin64/usr/local/bioinfo.json" ;;
#      cygwin64)   latest_package_url="${base_url}/cygwin64/cygwin64.json"          ;;
#    esac
#  fi

  if [[ -z ${latest_package_url:+set} ]]; then    # defaults needed for old releases without version and release json
    release_track="Beta"
    base_url="https://raw.githubusercontent.com/WGSExtract/WGSExtract-Dev/master/"
    latest_package_url="https://raw.githubusercontent.com/WGSExtract/WGSExtract-Dev/master/latest-release-Beta.json"
  fi

  # Retrieve a new (package) json file if the last one retrieved is not appropriate or existing
  if [[ ! -e latest.json || "$latest_package_url" != "$saved_package_url" ]]; then
    # if [[ ${release_track} == "Dev" ]]; then
    #  echo "Retrieving the WGS Extract latest available release file (track ${release_track}, package $1)."
    # else
      echo "Retrieving the WGS Extract latest available release file (track ${release_track})."
    # fi
    ${curlx} -o latest.json "$latest_package_url"
  fi

  export release_track base_url latest_package_url     # Provides source definition for latest.json file content
  # echo "get_latest_json: track $release_track, URL $latest_package_url"
}
export get_latest_json

declare latestVer latestDate latestURL
get_latest_release_info () {      # $1 is the package name (key in json)
  get_latest_json "$1"     # Returns release_track, latest_package_url and latest.json file for a given package

  latestDate=""
  if [ -e latest.json ]; then
    latestVer=$(jq -r .\""$1"\".\"version\" latest.json)
    latestDate=$(jq -r .\""$1"\".\"date\" latest.json)
    latestURL=$(jq -r .\""$1"\".\"URL\" latest.json)        # Actual pointer to this package ZIP archive for version
  fi

  if [[ -z ${latestDate:+set} ]]; then      # jq above may return "" if error in file content and format
    latestVer=-1
    latestDate="unk"
    latestURL="URL_unknown"
    echo "*** ERROR: No latest version info for package $1 in track ${release_track}"
  else
    echo "Found package $1 latest available version ${latestVer}, date ${latestDate} for track ${release_track}"
  fi

  export latestVer latestDate latestURL
  # echo "get_latest_release: ver $latestVer, date $latestDate, URL $latestURL"
}
export -f get_latest_release_info

declare currentVer currentDate currentURL
read_current_release_info () {      # $1 is the package name (key in json), $2 is the json file

  currentDate=""
  if [ -e "$2" ]; then    # May assign "" if error in file content or format
    currentVer=$(jq -r .\""$1"\".\"version\" "$2")
    currentDate=$(jq -r .\""$1"\".\"date\" "$2")
    currentURL=$(jq -r .\""$1"\".\"URL\" "$2")
  fi

  if [[ -z ${currentDate:+set} ]]; then
    currentVer=0
    currentDate="unk"
    currentURL="URL_unknown"
    echo "*** ERROR: No current version info for package $1 in local file $2"
  else
    echo "Found package $1 current installed version ${currentVer}, date ${currentDate}"
  fi

  export currentVer currentDate currentURL
  # echo "read_current_release: ver $currentVer, date $currentDate, URL $currentURL"
}
export -f read_current_release_info

# Special version for old Alpha 4m / 4.33 release that introduced version.json files but diff file convention
read_current_release_info_4m () {      # $1 is the package name (key in json), $2 is the json file
  currentDate=""
  pname="$1"
  if [[ $1 == "program" ]]; then
    pname="wgse"
  elif [[ $1 == "tools" ]]; then
    pname="localtools"
  fi
  if [ -e "$2" ]; then    # May assign "" if error in file content or format
    currentVer=$(jq -r .\""$pname".version\" "$2")
    currentDate=$(jq -r .\""$pname".date\" "$2")
    currentURL=$(jq -r .\""$pname".URL\" "$2")
  fi

  if [[ -z ${currentDate:+set} ]]; then
    currentVer=0
    currentDate="unk"
    currentURL="URL_unknown"
    echo "*** ERROR: No current version info for package $1 in local file $2 (Alpha 4m)"
  else
    echo "Found package $1 current installed version ${currentVer}, date ${currentDate} (Alpha 4m)"
  fi
  # echo "read_current_release_4m: ver $currentVer, date $currentDate, URL $currentURL"
}
export -f read_current_release_info_4m

declare newTrack
check_release_track() { # $1 is the new release.json file, existing is found in ./release.json
  newTrack=true     # Default is to force a new release.json file to be installed
  if [ -e release.json ] && [ -e "$1" ]; then
    currentTrack=$(jq -r .\"release\".\"track\" release.json)
    newTrack=$(jq -r .\"release\".\"track\" "$1")
    if [ "$currentTrack" == "$newTrack" ]; then
      newTrack=false
    fi
  fi
  export newTrack
}
export -f check_release_track

# ------------------------------------------------------------------------------------------------
# Install / Upgrade WGSE package system
#  Could potentially be local to zinstall_common.sh but zinstall_stage2windows.sh needs also

declare replace upgrade
install_or_upgrade() {    # $1 is package
  zipfile="$1.zip"        # Temporary zip file name after download amd before extract
  destpath="."            # Where to put the archive content

  # In install directory so paths relative to that
  case $1 in
    installer)
      verdir=scripts/               ;  longname="Installer"  ;;
    program)
      verdir=program/               ;  longname="Program"   ;;
    reflib)
      find_reflibdir
      verdir="${reflibdir}"         ;  longname="Reference Library" ;;
    tools)
      verdir=jartools/              ;  longname="Local Tools"  ;;
    bioinfo)    # Windows stage2 installer only; path local to the cygwin environment
      verdir=cygwin64/usr/local/    ;  longname="Bioinformatic Tools"
      destpath=$(cygpath -u cygwin64/usr)  ; zipdir="local"  ;;   # Override defaults
    cygwin64)   # Windows only (not really handled here; See install_windows.bat)
      verdir=cygwin64/              ;  longname="Cygwin64"  ;;
    *)
      echo "ERROR -- unknown package type $1"
      return ;;
  esac

  latestVer=0  ;  currentVer=0      # Zero out for safety from previous runs
  # Get latest available version info
  get_latest_release_info "$1"    # Sets three latest* variables & release_trackn
  if (( latestVer <= 0 )); then
    echo "*** WGS Extract package ""${longname}"" latest available version not found; no update possible."
    export replace=false  ;  export upgrade=false      # defaults going in; do / did nothing
    return
  fi
  # So latestVer valid at 1 or greater; so latestURL must be good and pointing to a valid package

  # Compare current to latest release info; decide what needs to be done.
  replace=false  ;  upgrade=false      # defaults going in; do / did nothing
  if [[ -d "$verdir" ]]; then
    current_json="${verdir}$1.json"         # Originally named version.json, now $package.json
    alpha4m_json="${verdir}version.json"    # Special for Alpha 4m only that has version.json files

    # echo current_json "$current_json" , alpha4m_json "$alpha4m_json"
    if [[ -e "$current_json" ]]; then
      read_current_release_info    "$1" "$current_json"      # Sets three current* variables; new package.json file
    elif [[ -e "$alpha4m_json" ]]; then
      read_current_release_info_4m "$1" "$alpha4m_json"      # Sets three current* variables; older version.json file
    else
      # echo "No JSON file found for package $1. Setting version to 0"
      currentVer=0    # Must be previous WGS Extract v3, v2, etc with no json file; so force upgrade

    fi
    if ((currentVer != 0)) ; then
      vmesg="v${currentVer}"
    else
      vmesg="older release"
    fi
    # Note: when copy latest installer over old release; adds scripts/installer.json. So installer will not upgrade.

    # Special exception.  We changed reflib version because we pulled the scripts out and put them into program pack.
    # So version reflects actual blob content of reference/ folder.  Version 35 is now known as version 5.
    if [ "$1" == "reflib" ] && ((33 <= currentVer && currentVer <= 35 && latestVer < 10 )); then
      checkVer=0            # Force upgrade to new, lower version number that is later
    elif [ "$1" == "reflib" ] && ((33 <= latestVer && latestVer <= 35 && currentVer < 10 )); then
      checkVer=latestVer    # We renumbered back to 5 from 35; so leave alone as <10 is newer than 33-35
    else
      checkVer=$currentVer  # Allows set and check of 0 override without wiping out currentVer value for reporting
    fi

    if (( checkVer < latestVer )); then
      start_mesg="*** UPGRADING WGS Extract ""${longname}"" from ${vmesg} ..."
      end_mesg="... finished upgrading WGS Extract ""${longname}"" to v${latestVer}."
      replace=true
      upgrade=true
    else
      echo "*** WGS Extract ""${longname}"" v${currentVer} is already installed and the latest available."
    fi
  else
    # We assume that if no json file then no "upgrade"; only new / "replace" (Installer will be latest)
    start_mesg="*** INSTALLING WGS Extract ""${longname}"" v${latestVer} ..."
    end_mesg="... finished installing WGS Extract ""${longname}"" v${latestVer}."
    replace=true
  fi

  # If local package needs to be replaced, do it now.
  if "$replace" ; then
    echo "$start_mesg"

    ${curlx} -o "$zipfile" "$latestURL"     # Get the package from the internet to a local .zip

    if [[ -e $zipfile ]]; then              # With downloaded package ...

      case $OSTYPE in                       # OS specific uncompress
        msys* | cygwin*)
          powershell Expand-Archive -LiteralPath "$zipfile" -DestinationPath "$destpath" -Force  ;;
        darwin*)    7zz x -tzip -y "$zipfile" >/dev/null  ;;
        linux*)     7z  x -tzip -y "$zipfile" >/dev/null  ;;
      esac
      if [ ! -d ${destpath}/${zipdir} ]; then
        echo "*** ERROR: ${zipdir} directory not created when expanding ${zipfile} during package ""$1"" install"
        return
      fi
      rm -f "$zipfile"

      case $1 in    # Handle package specific installation steps with unzipped files
        installer)
          # We generally will only replace the local release.json file if we change tracks (or it does not exist yet)
          check_release_track "$zipdir/release.json"    # Sets global $newTrack if changing release track
          if [[ ! -e release.json || $newTrack ]]; then
            if [[ -e release.json ]]; then
              mv release.json release_oldtrack.json     # Preserve old file in case dev changed to use their own URLs
            fi
            mv "${zipdir}/release.json" .
          fi

          # To force an update release.json even if the same track, anytime release-override.json stored in the archive
          if [ -e "${zipdir}/release-override.json" ]; then
            if [[ -e release.json ]]; then
              mv release.json release_overridden.json   # Preserve old file in case dev changed to use their own URLs
            fi
            mv "${zipdir}/release-override.json" release.json
          fi

          { # *** Critical code area ; this script changes after cp. So force a pre-read of critical section
          \cp -rf ${zipdir}/* .         # Simply copy everything including this script (overwrite)
          chmod a+x ./*.bat ./*.sh ./*.command scripts/*.sh   # Should be set in zip; but just in case

          # Ending this routine here as we are about to restart in a new script
          rm -rf ${zipdir} || true
          echo "$end_mesg"

          echo
          echo '****************************************************************************************************'
          echo '*** Need to restart the WGS Extract v4 installer with the upgraded installer scripts.'
          echo '****************************************************************************************************'
          echo
          read -r -p 'Press enter / return to restart the WGSE installer ...'
          case $OSTYPE in
            msys* | cygwin*)  $(cygpath "$COMSPEC") /c Install_windows.bat  ;;
            darwin*)          ${bashx} Install_macos.command  ;;
            linux*)           ${bashx} Install_ubuntu.sh  ;;
          esac
          exit 10     # Exit this script with status 10 as we already restarted the new one
          } ;;   # *** End of critical code area

        program)
          \cp -rf ${zipdir}/* .      # Includes program/, some scripts/, open_source_licenses/, WGSExtract*, Library*
          chmod a+x ./*.bat ./*.sh ./*.command scripts/*.sh  ;;

        reflib)
          [ ! -d "${reflibdir}" ] && mkdir "${reflibdir}"   # Used to just copy zipdir which created reference/
          pushd ${zipdir}/reference > /dev/null ;  \cp -rf ./* "${reflibdir}"  ;  popd > /dev/null

          \cp -f "${reflibdir}00README_genomes.txt" "${reflibdir}genomes/"  # Install new readme in genomes/ subdir
          \rm -f "${reflibdir}genomes/genomes.csv" || true   ;;  # Want to force new seed_genomes.csv

        tools)
          \cp -rf ${zipdir}/* .     # Simply copy everything; allows easy additions without modifying here
          chmod a+x FastQC/*.bat FastQC/fastqc  ;;

        bioinfo)
          ;;        # Expands (and removes) zip file directly into cygwin64/usr/local so nothing to do

        cygwin64|*)
          echo "*** ERROR - Unexpected package type $1 during Install or Upgrade"  ;;
          # cygwin64 is handled in the Install_windows.bat file before it brings in the BASH environment
      esac

      rm -rf ${zipdir} || true
      echo "$end_mesg"
    else
      echo "*** FAILURE when trying to download the package ""${longname}"" v${latestVer}"
    fi
  fi

  export replace upgrade
#  echo "install or upgrade $1 $2: replace $replace, upgrade $upgrade"
}
export -f install_or_upgrade

# ------------------------------------------------------------------------------------------------
# Reference Genome Library Seed file load
#  Code needed in get_and_process_refgenome.sh and zlibrary_common.sh so placed here
#  Seed genomes.csv has one row per genome; columns are:
#   Python Genome Code, Final File Name, Downloaded File Name, URL, Library command menu string, SN Cnt, SN Name, Descr
# Note: very order dependent here and in the genomes.csv file.  Also, true CSV with each field in double quotes and
#  comma separated (no commas or double quotes allowed in fields).

declare -a pytcode finalf initf gurl menopt sncnt snnam descr
read_genomes_file () {
  # Start with seed file from release if non-existent
  [ ! -f "${reflibdir}genomes/genomes.csv" ] && cp -f "${reflibdir}seed_genomes.csv" "${reflibdir}genomes/genomes.csv"

  # Read in Genomes.csv and transpose into array of column strings
  declare -a row cols
  while IFS="," ; read -ra row ; do
    for ((i=0 ; i < ${#row[@]}; i++)) ; do
      cols[$i]+="${row[$i]},"     # Rebuild rows by appending to column array strings; reinsert comma separator
    done
  done < "${reflibdir}genomes/genomes.csv"

  # This order (of array index / columns) must match the order in the genomes.csv file; first row is header / titles
  IFS=","  read -ra pytcode <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[0]}")   # Python Ref Gen Code eg hs37
  IFS=","  read -ra source  <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[1]}")   # Server source of file (NIH, EBI, etc)
  IFS=","  read -ra finalf  <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[2]}")   # Final reference file name
  IFS=","  read -ra initf   <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[3]}")   # Initial, downloaded file name
  IFS=","  read -ra gurl    <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[4]}")   # URL to download ref genome
  IFS=","  read -ra menopt  <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[5]}")   # BASH Library menu option
  IFS=","  read -ra sncnt   <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[6]}")   # Seq Name Count (# SNs)
  IFS=","  read -ra snnam   <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[7]}")   # Seq Name Format (Chr, Num, Acc)
  IFS=","  read -ra descr   <<< $(sed -e 's/"//g' -e 's/,$//' <<<"${cols[8]}")   # Description
  export pytcode source finalf initf gurl menopt sncnt snnam descr
}
export -f read_genomes_file

filebad() {  # file name, min size to check
  # Check if file bad: does not exist or is smaller than the minimum size allowed
  if [ -f "$1" ] ; then
    case $OSTYPE in
      msys* | cygwin* | linux*)  size=$(stat -c %s "$1") ;;
      darwin*)                   size=$(stat -f %z "$1") ;;
    esac
    (( size < $2 ))
  else
    true
  fi
}

get_and_process_refgenome() {   # index to genomes.csv arrays
  echo
  echo "Downloading and Processing ${descr[$1]}"
  cd "${reflibdir}"genomes/

  [ -f "${finalf[$1]}" ] && rm -f "${finalf[$1]}"
  [ -f "${initf[$1]}" ] && rm -f "${initf[$1]}"

  # To get around exclanatiobs embedded in MS Onedrive URLs; no spaces in params so no escaped double quotes needed
  IFS=" " read -r -a cmd <<< "${curlx} -o ${initf[$1]} ${gurl[$1]}"
  "${cmd[@]}"

  if filebad "${initf[$1]}" 500000000 ; then
    echo "*** Error downloading ${initf[$1]}"
    [ -f "${initf[$1]}" ] && rm -f "${initf[$1]}"
    return
  fi

  ${bashx} "$process_refgenomes" "${initf[$1]}"

  if filebad "${finalf[$1]}" 500000000 ; then
    echo "*** Error processing ${finalf[$1]}"
    [ -f "${finalf[$1]}" ] && rm -f "${finalf[$1]}"
    [ -f "${initf[$1]}" ] && rm -f "${initf[$1]}"
    return
  fi

  echo "Success installing ${pytcode[$1]}, ${finalf[$1]}, SN Cnt: ${sncnt[$1]}, SN Style: ${snnam[$1]}"
}
export -f get_and_process_refgenome
