#!/usr/bin/env bash
#
# Superset of function in the old v3 get_and_process_refgenomes.sh and UI of the v4 zlibrary_common.sh Library
#  command and get_and_process_refgenome.sh file. Simplified handling by relying on the genomes.csv file in the
#  reference/genomes folder. Due to returning to this external single file mechanism of v3, added a command line mode.
#  Expanded functionality to let developers verify various aspects and subsets of function; including using local files
#  already downloaded.  See usage below for more details.
#
# Called from the Library* user command and python referencelibrary.py. Runs standalone (direct BASH invocation).
#
# Part of the main program package in WGS Extract (https://wgsextract.github.io/) to support the Ref Genome Library
# Copyright (c) 2021-23 Randolph Harr
#

declare local_dir="/cygdrive/d/dna/reference_genomes_orig"    # A default path for the source of local file(s)
# declare -ir skip=6               # Skip first 6 genomes.csv entries that dup next 6 (NIH vs EBI) for final files
declare -ir adjust=5             # Number of Library Mode menu items before regular list begins (set in library mode)

read -r LINES COLUMNS < <(/bin/stty size)   # Set width of screen (characters) for table "list" (used by column cmd)

# Note: If (1) in keep mode, (2) initf and finalf file names are the same, (3) final file exists, and (4) the final
#       fileis smaller than the downloaded one; then curl thinks the file did not finish downloading and simply
#       downloads the remaining bytes. Thus corrupting the final file already there and having a corrupted initial file.

# Note: Have to have separate PREP call because if in stage-centric mode, then a single call to PROC handles all
#       files at once.  So need all files PREP'ped before calling that single PROC call. Otherwise, would have to
#       detect this in PROC stage call and handle the looping there accordingly.

usage() {
  echo
  printf "Usage: %s [ \"force\" ] { command } [ stages ] [ get_modifier ] \n" "${0##*/}"
  printf "  to GET, VERIFY, PREP_and_PROCESS and REVERIFY reference genomes in the WGS Extract ref library.\n"
  $debugmode && return
  echo
  printf "     \"force\"     Force delete a file before replacing.  Default is to keep verified files.\n"
  echo
  printf "  * A Required command from one of the following (no default):\n"
  printf "      file       a single genomes.csv file entry to handle (\"initf\" or \"finalf\" file name),  OR\n"
  printf "      index      a single genomes.csv file entry (numeric index; ignores NIH/EBI setting),   OR\n"
  printf "     \"all\"       handle ALL entries in the genomes.csv file in file-centric mode,            OR\n"
  printf "     \"stage\"     \"all\" but stage-centric mode (all files in a stage before the next stage),  OR\n"
  printf "     \"list\"      list ALL the installed reference genomes,                                   OR\n"
  printf "     \"verify\"    verify ALL the installed reference genomes (same as \"all pveronly\"),        OR\n"
  printf "     \"clean\"     delete ALL \"process\" created files from the reference library,              OR\n"
  printf "     \"clean_all\" \"clean\" ALL including retrieved reference genomes,                          OR\n"
  printf "     \"library\"   loop through a user menu (the Library.xxx command)\n"
  echo
  printf "  * An Optional limit to the stages processed (default is all stages):\n"
  printf "     \"get\"       Only retrieve and verify file(s); do not prep, process or reverify them,      OR\n"
  printf "     \"proc\"      Only prep, process and verify existing files; do not retrieve ,               OR\n"
  printf "     \"getonly\"   Only retrieve files,                                                          OR\n"
  printf "     \"gveronly\"  Only verify already retrieved files,                                          OR\n"
  printf "     \"proconly\"  Only prep-and-process already retrieved files,                                OR\n"
  printf "     \"pveronly\"  Only verify already processed files.\n"
  echo
  printf "  * An Optional modifier to GET stage handling:\n"
  printf "     \"NIH\"       URL download mode and use the NIH and not the EBI server, where possible,     OR\n"
  printf "     \"EBI\"       URL download mode and use the EBI and not the NIH server, where possible.     OR\n"
  printf "     \"local\"     copy file(s) from a \"local_dir\"; do not download,                             OR\n"
  printf "      dir        \"local\" copy mode using the absolute path \"dir\" as the \"local_dir\".\n"
  printf "    The default is URL download mode; no preferred server (both if all).\n"
  printf "    The default \"local_dir\" is %s\n" "$local_dir"
  echo
  printf "  * File-centric is all stages per file before the next file (default).\n"
  printf "    Stage-centric is all files per stage before the next stage.\n"
  printf "  * File and dir(ectory) names should be double-quoted to prevent shell globbing.\n"
  printf "  * The \"library\" command menu number is 5 more than than the genomes.csv \"index\" used here.\n"
  echo
}


#------------------------------------------ SETUP WGSE ENVIRONMENT --------------------------------------------------

# Find the installation directory (normally part of zcommon.sh but when run standalone, cannot find that script)
WGSEDIR=$(/usr/bin/dirname "${BASH_SOURCE[0]}")   # Get the script location to determine install directory
WGSEABS=$(\cd "$WGSEDIR" || true; \pwd -P)                # By cd'ing to it, resolve any aliases and symlinks
[[ $(/usr/bin/basename "${WGSEABS}") == "scripts" ]] && WGSEFIN=$(/usr/bin/dirname "${WGSEABS}") || WGSEFIN="${WGSEABS}"

# For local shellcheck, declare variables and functions exported from zcommon.sh that we later use
declare -f find_reflibdir
declare reflibdir

declare -f read_current_release_info
declare currentVer
declare currentDate
declare libver=0                    # Local to support banner reporting Reference Library version

declare -f read_genomes_file
declare -a source finalf initf gurl menopt descr md5f md5i sizef sizei  # sortord pytcode sncnt snnam build class
declare -i genlen=0                 # Local to support # entries in genomes.csv file (including header)

declare oWGSEFIN                    # OS Native form of isntalltion path
declare curlx bashx rmfx cpfx mvfx
declare process_refgenomes
declare -i success fail

if ! source "${WGSEFIN}/scripts/zcommon.sh" dummy ; then
  echo "*** ERROR sourcing scripts/zcommon.sh in ${BASH_SOURCE[0]}. Exiting."
  exit 1
fi


# ------------------------------------------ PROGRAM MODE SETTINGS --------------------------------------------------

# See zcommon.sh for our special note on pseudo-boolean variables in BASH as used here

declare library=false         # if true, library mode (user queries in UI); else command_line mode

declare stagecentric=false    # if true, all files for a stage handled before doing next stage; else file-centric

declare keepmode=true         # Reuse verified files if true (default). Otherwise, force delete files before a stage.
declare debugmode=false       # internal option (undocumented); turn on verbose debug messages in keep mode
declare debugstr=""           # string version of debug mode to pass to process_refgenomes.sh (only option is keep)

declare -a stages=(get gver prep proc pver)   # Five stages of possible processing on a file (override with parameters)
                              # Note: user only given option for 4. Prep always performed before Proc.

declare -i oneindex=0         # the numeric index of the file in the genomes.csv table to process
declare all=false             # if true, then all files (oneindex=0) (convenience over oneindex alone)

declare notsource="none"      # if NIH, then (not) EBI-Alt. If EBI, then (not) NIH-Alt (during get file only)
declare localmode=false       # if true, use local file(s). False means use the URL. (only applies to get / retrieve)
# declare local_dir             # set at top of this file ; changed by "dir" spec

declare localfile=""          # global path for processing a local oneindex entry (set as needed in localmode)

declare -a notverinit=()      # To avoid running md5sum multiple times; If false, then md5 checked and good. Test of an
declare -a notverifinal=()    # unset variable returns true; so inverse logic. Need array for stage-centric processing.

declare -a prepf=()           # Mod to initf[index] after prep so used in process call to generate finalf[index]


#-------------------------------------------- LOCAL UTILITY SUPPORT ------------------------------------------------

is_index_good() {
  # Check $1 index into genomes.csv array for valid range; report to stderr with $2 routine name if bad
  local mesg=""

  # shellcheck disable=SC2236
  if [[ -n $1 ]] && (( $1 > 0 && $1 < genlen )) ; then
    return $success        # Index is good
  else
    [[ -n $2 ]] && mesg=" inside $2" || mesg=""
    echo "*** Error: bad index \"$1\" to genomes.csv array${mesg}." >&2
    return $fail       # Index is bad; reported as such
  fi
}


is_file_good() {
  # $1 is index, $2 init mode? (vs final), $3 report error?
  # File good IF: verified already, it exists, is correct size, the expected md5sum matches and (if final) processed
  # SIDE EFFECT: Sets global array notverxxx[] to false for index IF file is good
  local verif veriff verifi md5chk report init final md5so base finalwgse notverified filefn prindex mentry
  declare -i sizechk sizescaled sizechkscaled lindex actual_size init_size final_size

  $debugmode && echo -n "*** DEBUG: is_file_good: index: $1, init? $2, report? $3"

  # For clarity; assign the parameters to a descriptive variable
  declare -ir index=$1
  $2 && init=true   || init=false
  $2 && final=false || final=true
  $3 && report=true || report=false

  # Should not be the case but check and report anyway
  if ! is_index_good $index "is_file_good" ; then
    $debugmode && wcho "" && echo "*** ERROR: internal error due to call with bad index (is_file_good)"
    return $fail
  fi

  # Create local variables based on index and mode (init or final)
  $final && md5chk=${md5f[$index]} || md5chk=${md5i[$index]}
  $final && sizechk=${sizef[$index]} || sizechk=${sizei[$index]}
  $final && filefn="${finalf[$index]}" || filefn="${initf[$index]}"
  $final && notverified="${notverifinal[$index]}" || notverified="${notverinit[$index]}"
  $notverified && verif=false || verif=true   # testing unset array variable is true so use not notverified
  ${notverifinal[$index]} && veriff=false || veriff=true
  ${notverinit[$index]}   && verifi=false || verifi=true

  $debugmode && echo ", verified? $verif, file: $filefn"

  $final && mesg="processed file" || mesg="retrieved file"

  # We adjust the printed index for (a) library mode menu and (b) in all and final mode,
  #   use both NIH & EBI index for printout because we skip EBI entries so we only check unique final values
  $library && lindex=$(( index + adjust )) || lindex=$index
  prindex="$lindex)"
  if $all && $final ; then
    if [[ ${source[$index]} == "EBI-Alt" ]] ; then
      $debugmode && echo "*** INTERNAL ERROR: Call to verify in all & final with EBI-ALT set"
      return $fail
    elif [[ ${source[$index]} == "NIH-Alt" ]] ; then
      # Need to find EBI-Alt and use both indices for output
      for i in "${source[@]}" ; do
        if [[ "${pytcode[$index]}" == "${pytcode[$i]}" && ${source[$i]} == "EBI-Alt" ]] ; then
          prindex="$lindex/$i)"
          break
        fi
      done
    fi
  fi

  # Simplify library mode menu entry to become the printed message reference model
  mentry=$(echo "${menopt[$index]}" | sed 's/ (NIH)//;s/ (EBI)//;s/ (\*)//;s/ 3x//')   # strip unneeded trailers

  # Check if file already verified; return if so (Debug mode reported value earlier)
  case $init$verifi$veriff in
    truetrue* | false*true)
      $report && echo "*** SUCCESS: Verified $prindex $mentry ($mesg)"
      return $success  ;;
    truefalsetrue)
      $report && echo "*** SUCCESS: Verified $prindex $mentry (processed file; not $mesg)"
      return $success  ;;
    falsetruefalse)
      $report && echo "*** ERROR: Unverified $prindex $mentry (found retrieved file; not $mesg)"
      return $fail  ;;
    *falsefalse)
      ;;              # Simply continue as not already verified (not meeded; left for clarity)
  esac

  # Check if file exists
  if [ ! -f "$filefn" ] && [[ ! $all || $debugmode ]] ; then
    # Shortcut to not report file missing when doing whole library
    $report && echo "*** ERROR: Unverified $prindex $mentry (missing $mesg)"
    return $fail

  fi

  # Check if file is the correct size
  case $OSTYPE in
    msys* | cygwin* | linux*)  actual_size=$(stat -c %s "$filefn") ;;
    darwin*)                   actual_size=$(stat -f %z "$filefn") ;;
  esac
  
  init_size="${sizei[$index]}"
  final_size="${sizef[$index]}"
  
  if (( actual_size != sizechk )) ; then

    if $final && (( actual_size == init_size )) ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (size, only retrieved file found)"

    # Will set verified if found final file in get routine earlier; so only get here if doing gver in stage mode
    elif ! $final &&  (( actual_size == final_size )) ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (size, only processed file found)"

    elif ! $final && (( actual_size < init_size )) ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (size, retrieved file partially downloaded?)"
      $report && echo "           Try downloading again to finish."

    else
      sizechkscaled=$((sizechk/1000000))
      sizescaled=$((actual_size/1000000))
      $report && echo "*** ERROR: Unverified $prindex $mentry (size mismatch, $mesg): expected $sizechkscaled MB, found $sizescaled MB."

    fi
    return $fail

  fi

  # Check md5sum of file (note: likely 3-5 second operation)
  md5so=$(md5sum "$filefn" | cut -d' ' -f1 | tr -d '\n')  #  cut is passing the trailing \n added by md5sum?
  if [[ "$md5so" != "$md5chk" ]] ; then

    # Special case where final file size is correct but md5sum is of an unprocessed, downloaded initial file.
    if $final && [[ "$md5so" == "${md5i[$index]}" ]] ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (md5, only retrieved file found)"

    elif ! $final && [[ "$md5so" == "${md5f[$index]}" ]] ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (md5, only processed file found)"

    else
      $report && echo "*** ERROR: Unverified $prindex $mentry (mD5sum mismatch, $mesg):" &&
                 echo "           expecting \"$md5chk\"" &&
                 echo "           found     \"$md5so\""

    fi
    return $fail

  fi

  # If checking final file, check that it was properly processed (last step is to create the .wgse catalog file)
  if $final ; then
    base=$(echo "$filefn" | sed "s/.fasta.gz//;s/.fna.gz//;s/.fa.gz//")
    finalwgse="${base}.wgse"
    if [ ! -f "$finalwgse" ] || [ "$filefn" -nt "$finalwgse" ] ; then
      $report && echo "*** ERROR: Unverified $prindex $mentry (bad catalog entry)"
      return $fail

    fi
  fi

  # Todo should we add an hstfile call to (a) verify is a FASTA and (b) if final, verify is bgzip ?

  $report && echo "*** SUCCESS: Verified $prindex $mentry ($mesg)"
  $final && notverifinal[index]=false || notverinit[index]=false     # Meaning, now verified
  return $success
}


prep_filename() {
  # $1 is index into the genomes.csv file; output is prepfn[$1] name in global scope
  # prepfn result may need to be different from genomes.csv initfn specification because:
  #   (a) process_refgenomes.sh may change the file name (e.g. add .gz after compression, change .zip to .gz. etc); or
  #   (b) for clarity, we set the initf[] to the cloud server common name but use a more succint, recognizable
  #       final file name in the reference/genomes directory (that is unique)
  # Using initfn, create necessary global prepfn for process_refgenomes so correct finalfn is created
  # Errors are considered internal errors due to improper genomes.csv entries and printed to stderr
  # todo why using global scope initfn and finalfn; why not genomes.csv array entries based on index?
  local bad finalext final2ex finalbas initext init2ex initbas finalfn initfn
  declare -ir index=$1

  # For clarity and simplicity below
  finalfn="${finalf[$index]}"   ;   initfn="${initf[$index]}"

  $debugmode && echo "*** DEBUG: prep_filename: index $index initfn $initfn, finalfn $finalfn"

  # First, a quick sanity check of the finalfn name for a valid extension of .[fa|fna|fasta].gz
  finalext="${finalfn##*.}"    # finalfn last extension (must be .gz)
  finalbas="${finalfn%.*}"     # finalfn name less last extension
  final2ex="${finalbas##*.}"   # finalfn 2nd to last extension (must be [.fa|.fna|.fasta])
  
  bad=false
  if [[ "$finalext" != "gz" ]]; then
    bad=true
  fi
  case ${final2ex} in
    fa|fasta|fna)             ;;    # finalfn properly ends in  .[fa|fna|fasta].gz ; do nothing
    *)            bad=true    ;;
  esac
  if $bad ; then
    echo "*** ERROR: finalf from genomes.csv must end in .[fa|fna|fasta].gz  (${finalfn})" >&2
    return $fail
  fi

  prepf[index]=$initfn    # Initial default; same as initfn

  if [[ "$initfn" == "$finalfn" || ! -f "$initfn" ]]; then
    # initfn and finalfn filenames are the same so we know finalfn is properly named. Nothing to do.
    #  OR initfn does not exist (just return)
    $debugmode && echo "*** DEBUG: prep_filename: initfn $initfn does not exist or is the same name as finalfn"
    return $success
  fi

  # Change the initfn name (if at all) so we end up with the finalfn after the process_refgenomes() call
  # initfn exists and must end in .fa|fna|fasta|zip|7z|bz|bz2|gz or is in error
  $debugmode && echo "*** DEBUG: prep_filename: prep initfn $initfn for use in processing"
  initext="${initfn##*.}"    # initfn last extension (must be .gz)
  initbas="${initfn%.*}"     # initfn name less last extension
  init2ex="${initbas##*.}"   # initfn 2nd to last extension (must be [.fa|.fna|.fasta])

  # Key below is to set prepf[] entry AND rename the file to the prepf[] from initfn
  case ${initext} in

    # Both files end in .gz but are different base names; rename the downloaded initfn file to finalfn
    gz)
      prepf[index]="$finalfn"
      $mvfx "$initfn" "${prepf[index]}"
      ;;

    # initfn is compressed but not in .gz format; process_refgenome will change compression and extension to .gz
    zip|7z|bz|bz2)
      # initfn base name != finalfn base name so need initfn to be finalfn base with initfn ext
      if [[ "${initbas}" != "${finalbas}" ]] ; then
        prepf[index]="${finalbas}.${initext}"
        $mvfx "$initfn" "${prepf[index]}"

      # else both filenames less extension are identical and initfn ends in non .gz compression extension
        #   process_refgenome will take care of this name and format change
      fi
      ;;

    # initfn is mot compressed; process_refgenome() will compress and add a .gz extension. So initfn+.gz == finalfn
    fa|fasta|fna)
      if [[ "${initfn}.gz" != "${finalfn}" ]]; then
        prepf[index]="${finalbas}"
        $mvfx "$initfn" "${prepf[index]}"      # note: handles if init2ex differnt than final2ex

      # else initfn ends in uncompressed extension and with .gz added is identical to finalfn
        #   process_refgenome will take care of this name and format change
      fi
      ;;

    # Bad extension on initfn
    *)
      echo "*** ERROR: initf extension (${initfn##*.}) must end in [fa|fasta|fna|zip|7z|bz|bz2|gz]" >&2
      return $fail
      ;;
  esac

  # return prepfn variable array element in caller scope; needed for single file or file-centric processing
  $debugmode && echo "*** DEBUG: prep_filename (finish): initfn \"$initfn\", prepfn[$index] \"${prepf[$index]}\", finalfn \"$finalfn\""

  return $success
}


handle_stage() {
  # Handles just one index ($1) in one of five stages ($2).  Called from a 2 nested loop genome handler
  # When index is zero, $all must be set and stage ($2) must be "proc"
  # Uses fixed stage names, genomes.csv arrays, mode settings.
  local report init final mentry stage dirmode finalfn initfn md5fd md5id prepfn action cmd
  declare -i lindex

  # For clarity and simplicity below
  declare -ir index=$1
  stage=$2

  init=true    ;  final=false
  report=true  ;  noreport=false

  { $all && $stagecentric ; } && dirmode=true || dirmode=false
  $library && lindex=$((index + adjust)) || lindex=$index

  #  Should never be an issue but just in case
  if ! is_index_good $index "handle_stage" ; then
    return $fail
  fi

  # Only set if valid index (just in case called with zero; should not be but ...)
  finalfn="${finalf[$index]}"   ;   initfn="${initf[$index]}"
  md5fd="${md5f[$index]}"       ;   md5id="${md5i[$index]}"
  prepfn=${prepf[$index]}       #   notverified() left alone
  mentry=$(echo "${menopt[$index]}" | sed 's/ (NIH)//;s/ (EBI)//;s/ (\*)//;s/ 3x//')   # strip redundant trailers

  if $debugmode ; then
    { $dirmode && [ "$stage" == "proc" ] ; } && echo "*** DEBUG: handle_stage: Process (dir mode)" || \
                                                printf "*** DEBUG: handle_stage: %-5s  %d\n" "$stage" "$index"
  fi

  # ! $dirmode && report=true || report=false     # Do not report if stagecentric and all mode
  $debugmode && report=true                     # Always report in Debugmode

  case $stage in

    #------------------------------------------------------------------------------------------------------------
    get)
      # Sometimes the curl fails prematurely.  if keepmode on, then let curl pickup where it stopped and finish
      # when run again.  BUT if already finished downloading (exact expected size) and keepmode on,
      # then rely on the md5 check to determine if we can reuse the file.  BUT, when the finalf exists and correct,
      # then skip as already downloaded, processed, and verified.

      if $keepmode ; then

        # If final correct file exists; reuse as if initf good if in keepmode
        if is_file_good "$index" "$final" "$noreport" ; then   # Do not report this check
          # notverifinal[$index]=false set in is_file_good call
          notverinit[index]=false  # Also set verified for init as reusing final
          $report && echo "*** SUCCESS: Reusing processed file $lindex) $mentry (get)"
          return $success

        # If init correct file exists; reuse as good if in keepmode
        elif is_file_good "$index" "$init" "$noreport" ; then
          # notverinit[$index] set false in is_file_good call
          $report && echo "*** SUCCESS: Reusing retrieved file $lindex) $mentry (get)"
          return $success

        # else  simply fall through and see if we can run the curl command again to restart and finish the download
        fi
      
      elif [ -f "$initfn" ] ; then    # and not in $keepmode ; may clear out a final file if the same name
        $rmfx "$initfn" || true    # Clear out any previous version before getting new one (! $keepmode)
        notverinit[index]=true    # Cleared file so clear verification

      fi

      $localmode && action="Copying local" || action="Downloading"

      # if localmode, then set localfile, verify it exists, and then perform the copy
      if $localmode ; then

        localfile="${local_dir}/${initfn}"   # Create full path to find file (use as a parameter later)
        if [ -f "$localfile" ] ; then
          if ! $cpfx "$localfile" "${reflibdir}genomes/${initfn}" ; then
            echo "*** ERROR: $action file failed: $lindex) $mentry"
            return $fail

          fi

        elif ! $all ; then
          echo "*** ERROR: Local file required but not found: $localfile"
          return $fail

        # When processing a whole directory (ALL), continue silently if a file is missing else ...
        fi

      #  URL Mode (not $localmode)
      else

        # Main curl command to download a reference genome; no spaces in params so escaped quotes not needed
        # Need read to form array to avoid shell interpretation of exclamations in MS Onedrive URLs
        IFS=" " read -r -a cmd <<< "${curlx} -o $initfn ${gurl[$index]}"
        $debugmode && echo "*** DEBUG: get: " "${cmd[@]}"
        if ! "${cmd[@]}" ; then  # Actual execute of curl command; curl has retry built in
          echo "*** ERROR: $action file failed: $lindex) $mentry"
          # Leave file, if it exists, as maybe another curl run can pick up where it left off.
          return $fail

        fi

      fi
      echo "*** SUCCESS: $action file: $lindex) $mentry"      # Always report success for an error free get
      ;;

    #------------------------------------------------------------------------------------------------------------
    gver)
      if ! is_file_good "$index" "$init" "$report" ; then
        # Reported error (or success) in is_file_good
        ! $keepmode && { notverinit[index]=true ; $rmfx "$initfn" || true ; }   # Clear out failed file
        return $fail

      fi
      ;;

    #------------------------------------------------------------------------------------------------------------
    prep)
      # Change initfn file so proc creates correct finalfn; returns success even if initfn does not exist
      # Have to split from proc) stage call as all and stagecentric means first proc call assumes all filenames prepped
      if ! prep_filename "$index" ; then     #  resets notverified
        $debugmode && echo "*** DEBUG: Prep_filename return $?"
        return $fail     # Internal error due to bad genomes.csv entries; already reported error inside function

      fi

      # Normally part of proc like is part of get. Checking if valid final file already exists.
      # We do it here because in dirmode, prep is called for each file whereas proc is called only once
      # finalfn may exist from previous run
      if $keepmode ; then
        if ! ${notverifinal[$index]} ; then
          # If already verified then simply return without message (likely reported in get call)
          return $success

        elif $all && [[ ${source[$index] == "EBI-Alt"} ]] ; then
          # If doing all, then ignore NIH-Alt as duplicated with EBI-Alt
          return $success

        elif is_file_good "$index" "$final" "$noreport" ; then    # even in debugmode, do not report error
          # Call to is_file_good will set notverifinal[$index] to False if checked out OK
          (( index <= skip*2 )) && prindex="$(( index - skip ))/${lindex}" || prindex="$lindex"
          $report && echo "***SUCCESS: Reusing verified processed file $lindex) $mentry (prep)"
          return $success

        fi

      elif ! $keepmode && [ -f "$finalfn" ] && [ "${prepf[$index]}" != "$finalfn" ] ; then
        ${bashx} "$process_refgenomes" "$debugstr" "clall" "$finalfn"   # Remove previous final files
        notverifinal[index]=true

      fi
      ;;

    #------------------------------------------------------------------------------------------------------------
    proc)
      # in directory mode process_refgenomes call (creates WGSE.csv, etc instead of simply appending to it)
      prepfn=${prepf[$index]}
      $debugmode && echo "*** DEBUG: proc: dirmode? $dirmode, notverified? ${notverifinal[$index]}, index $index, prepfn $prepfn"

      if $dirmode && (( index == 1 )) ; then    # Always do all if dirmode and first call; do not worry about what exists
        ${bashx} "$process_refgenomes" "$debugstr" "${reflibdir}genomes/"

       # Single file process_refgenones call (not needed if prepfn does not exist)
      elif ! $dirmode && [ -f "$prepfn" ]  && ${notverifinal[$index]} ; then
        ${bashx} "$process_refgenomes" "$debugstr" "$prepfn"    # NOOP when everything already run

      # else silently skip if all; report missing file only if single file mode during debug
      elif $debugmode ; then
        echo "*** DEBUG: Skipping $lindex) $mentry - file does not exist (prep)"

      fi
      return "$?"
      ;;

    #------------------------------------------------------------------------------------------------------------
    pver)
      if $all && [[ "${source[$index]}" == "EBI-Alt" ]] ; then   # Skip EBI-Alt final (identical to NIH-Alt) in ALL mode
        $debugmode && echo "*** DEBUG: ignore $index) $finalfn that we checked already (NIH vs EBI) (pver)"
        return $success

      elif ! is_file_good "$index" "$final" "$report" ; then
        ! $keepmode && ${bashx} "$process_refgenomes" "$debugstr" "clall" "$finalfn"
        return $fail

      fi
      ;;

  esac

}

#-------------------------------------------- MAIN PROCESSING CALL --------------------------------------------------
handle_refgenomes () {
  # $1 is index into genomes.csv entry to handle (ignored or 0 if $all is true)
  # Utilize all globals like keepmode, local_dir, localmode, stages(), etc
  local mesg first
  declare -i prindex outbeg outend inbeg inend error stagelen out in stage
  declare -i index=$1                             # For clarity below; not -r because will reuse inside loops

  # Setup initial banner message; list all stages set unless both get and proc; then simply state that
  if [[ " ${stages[*]} " =~ " get " && " ${stages[*]} " =~ " proc " ]] ; then
    mesg="Retrieve and Process"
  else
    mesg=""
    [[ " ${stages[*]} " =~ " get " ]] && { $localmode && mesg+="Copy" || mesg+="Download" ; }
    [[ " ${stages[*]} " =~ " gver " ]] && { [[ -z $mesg ]] && mesg+="Verify (retrieved)" || mesg+=", Verify (retrieved)" ; }
    [[ " ${stages[*]} " =~ " proc " ]] && { [[ -z $mesg ]] && mesg+="Process" || mesg+=", Process"; }
    [[ " ${stages[*]} " =~ " pver " ]] && { [[ -z $mesg ]] && mesg+="Verify (processed)" || mesg+=", Verify (processed)" ; }
    if [[ ! " ${mesg} " =~ ", "  ]] ; then
      first=${mesg:0:1}
      mesg="Only ${first,,[$first]}${mesg:1}"
    fi
  fi
  $library && prindex=$((index + adjust)) || prindex="$index"
  $all && mesg+=" Reference Genomes Library." || mesg+=" ${prindex}) ${descr[$1]}"

  echo
  echo "$mesg"

  if ! $all && ! is_index_good $index "handle_refgenomes" ; then
    return $fail
  fi

  # Lets clear out verified arrays first so forced to start from scratch like program startup
  notverinit=()
  notverifinal=()

  # mesg=$(echo "$mesg" | tr '[:upper:]' '[:lower:]' )    # Prep for later use in Success message

  # We call each individual stage per file via two nested loops.
  # We set the loop parameters depending on the mode we want
  # Outer loop cycles through file(s)  if file centric; else stage(s) if stage centric
  # Inner loop cycles through stage(s) if file centric; else file(s)  if stage centric
  # Single-file mode is run with a single index in file centric mode; all mode loops through all files in genomes.csv

  # Set loop start / stop limits depending on stage / file -centric and all / single file(s)
  stagelen=${#stages[@]}
  $debugmode && echo "*** DEBUG: all=$all, stagecentric=$stagecentric, genlen=$genlen, stagelen=$stagelen"
  # [[ $all && $stagecentric ]] && dirmode=true || dirmode=false
  case $all$stagecentric in
    truetrue)   outbeg=0      ;  outend=$(( stagelen  ))        # Stages (Stage Centric)  (dirmode)
                 inbeg=1      ;   inend=$(( genlen    ))  ;;    # Files
    truefalse)  outbeg=1      ;  outend=$(( genlen    ))        # Files (File Centric)
                 inbeg=0      ;   inend=$(( stagelen  ))  ;;    # Stages
    false*)     outbeg=$index ;  outend=$(( index + 1 ))        # File  (single file; not all; verified good index)
                 inbeg=0      ;   inend=$(( stagelen  ))  ;;    # Stages
  esac
  $debugmode && echo "*** DEBUG: outbeg=$outbeg, outend=$outend, inbeg=$inbeg, inend=$inend"

  error=0
  for (( out=outbeg; out<outend; out++ )) ; do
    for (( in=inbeg; in<inend; in++ )) ; do

      case $all$stagecentric in
        truetrue )            stage=$out ; index=$in  ;;    # dirmode
        truefalse | false* )  index=$out ; stage=$in  ;;    # file centric (all or single file)
      esac

      # Only break inner loop of stages if single file (not all mode) and stage error
      if ! handle_stage "$index" "${stages[$stage]}" && ! $stagecentric ; then
        # Error already reported in handle_stage
        error+=1
        break   # Break inner loop so no further processing on the file

      fi
    done

  done
  return $error
}


list_refgenomes() {
  local typind filelist unprocessed entry fentry mentry  finalfn initfni initfnj fentry mentry base finalwgse
  declare -i sec i j li lj

  noreport=false  # ; final=false  # Convenience in parameters to is_file_good below

  declare -ir hwidth=15 first=6   # adjust=5  - set globally at top
  # If file name longer than 2*hwidth+3 (3 dots), then drop middle of file name so get 2*hwidth len
  # When truncating file name, keep $first characters and last 2*hwidth-$first (note: no longer using file name)

  $library && typind="by menu entry" || typind="by genomes.csv index"

  echo
  echo -n "Listing installed Reference Genomes ($typind) "    # No newline as we print a . for each file found
  $debugmode && echo

  filelist=""
  unprocessed=false
  for ((i=1 ; i<genlen ; i++)); do   # genlen includes comment row 0 so <= not needed; skip EBI entries

    $debugmode && >&2 echo "*** DEBUG: in list: $i finalf ${finalf[$i]}"  # stdout redirected to column command ...

    (( i <= skip*2 )) && j=$(( i - skip )) || j="$i"  # Need to check NIH init file j that may be different than EBI i
    finalfn="${finalf[$i]}"
    initfni="${initf[$i]}"
    initfnj="${initf[$j]}"

    fentry=$(echo "$finalfn"      | sed 's/.fasta.gz//;s/.fna.gz//;s/.fa.gz//')      # strip known extension
    mentry=$(echo "${menopt[$i]}" | sed 's/ (NIH)//;s/ (EBI)//;s/ (\*)//;s/ 3x//')   # strip redundant trailers

    if [[ -f "$finalfn" || -f "$initfni" || -f "$initfnj" ]] ; then

      # Adjust string of file name so no more than 2*hwidth+3 characters wide ; always using (same) final file name
      if (( ${#fentry} > 2*hwidth+3 )) ; then
        sec=$(( -1 * (2*hwidth-first) ))         # last section = max length to be allowed minus skipped first
        fentry="${fentry:0:$first}...${fentry:$sec}"
      fi

      # adjust set globally at top of file; menu vs genomes.csv index. print indices for NIH and EBI files
      $library && li=$((i + adjust)) || li=$i
      $library && lj=$((j + adjust)) || lj=$j

      # When the initial files are different; then report both file indices
      (( i != j )) && prindex="${lj}/${li}" || prindex="$li"

      # Use library menu entry if library mode; else use final file name (possibly shortened)
      $library && entry="${prindex}) $mentry" || entry="${prindex}) $fentry"

      # if ! is_file_good "$i" "$final" "$noreport" ; then    # too time consuming; leave to verify command
      base=$(echo "$finalfn" | sed "s/.fasta.gz//;s/.fna.gz//;s/.fa.gz//")
      finalwgse="${base}.wgse"
      if [ ! -f "$finalwgse" ] || [ "$filefn" -nt "$finalwgse" ] ; then
        entry+=" ++"
        unprocessed=true
      fi

      # Create list with newlines that we pass to "column" at end
      filelist+="${entry}"$'\n'
      ! $debugmode && echo -n "."

    fi
  done
  echo      # Complete first line

  if [[ -n "$filelist" ]] ; then
    column <<< "$filelist"
    $unprocessed && echo "++ not cataloged"    # Only happens if filelist has content and one exists
  else
     echo "*** None Found"
  fi
}


# Save and restore is only needed for Library call; process args exits immediately after verifying the folder
verify_refgenomes() {
  local save_all save_stages save_stagecentric

  save_all="$all"                     ;  all=true
  save_stages=("${stages[@]}")        ;  stages=( pver )
  save_stagecentric="$stagecentric"   ;  stagecentric=false     # force file centric to get messages

  handle_refgenomes "0"       # When not all or not stage centric, pver prints individual success or error message

  all="$save_all"
  stages=("${save_stages[@]}")
  stagecentric="$save_stagecentric"
}


# ---------------------------------------- PROCESS SCRIPT PARAMETERS -----------------------------------------------
process_args() {
  # Pass in parameter list same as to script ("$@") ;  hassle is no shift available so must use our own loop
  # Make sure to call read_genomes_file before processing the args; so we can verify index and file name args

  # declare -a pargs=( "$@" )   # Convert args list to an indexable array; could use ${!i} and shift; but still need optN
  local opt="opt1"
  local valid=false
  local onefile=""
  local -i j

  # Loop through args (indicated by loop index i) ; mode for stage of arg processing is "opt"
  # while [[ $1 ]] will loop through elements of first arg ; scalar $1 in loop is that list element
  while [[ $1 ]] ; do

    # Could have let most fall through to the next case (;&) but simply use continue instead
    case $opt in

      ### optional [ test | debug ] (only)
      opt1)
        if [[ "$1" == "force" ]] ; then       # default is keep (not force) and reuse verified files
          shift ; keepmode=false

        elif [[ "$1" == "debug" ]] ; then     # Undocumented option; mutually exclusive with force
          shift ; debugmode=true  ;  keepmode=true  ;  debugstr="debug"   # Will report DEBUG mode is ON later

        fi

        opt="opt2" ; continue
        ;;

      ### Required main mode specification
      opt2)
        if [[ "$1" == "library" ]]; then
          library=true

        elif [[ "$1" == "clean" ]]; then
          ${bashx} "$process_refgenomes" "$debugstr" "clean"
          exit

        elif [[ "$1" == "clean_all" ]]; then
          ${bashx} "$process_refgenomes" "$debugstr" "clean_all"
          exit

        elif [[ "$1" == "list" ]]; then
          list_refgenomes
          exit

        elif [[ "$1" == "verify" ]]; then
          verify_refgenomes
          exit

        elif [[ "$1" == "all" ]]; then
          all=true  ;  stagecentric=false   ;  oneindex=0

        elif [[ "$1" == "stage" ]]; then
          all=true  ;  stagecentric=true    ;  oneindex=0

        # Check for valid index in genomes.csv; set onefile with the index
        elif [[ $1 =~ ^[0-9]+$ ]] ; then
          if ! is_index_good "$1" "process_args" ; then
            valid=false
            break     # Error reported in call above
          fi
          oneindex="$1"                   # Single file, good numeric index

        # Assume a single genomes.csv entry file name as the requied parameter (will not be blank)
        else
          onefile="${1##*/}"              # strip any pathname component; check validity later

        fi

        shift ; valid=true
        opt="opt3" ; continue
        ;;

      ### optional stage limits [ get | proc | getonly | gveronly | proconly | pveronly ] -- No prep on its own
      opt3)
        if [[ "$1" == "get" ]] ; then
          shift ; stages=( get gver )

        elif [[ "$1" == "proc" ]] ; then
          shift ; stages=( prep proc pver )

        elif [[ "$1" == "getonly" ]] ; then
          shift ; stages=( get )

        elif [[ "$1" == "gveronly" ]] ; then
          shift ; stages=( gver )

        elif [[ "$1" == "proconly" ]] ; then
          shift ; stages=( prep proc )

        elif [[ "$1" == "pveronly" ]] ; then
          shift ; stages=( pver )
        fi

        opt="opt4" ; continue
        ;;

      ### optional [ NIH | EBI | local | dir ] modifiers for "get" stage]  (last opts)
      opt4)
        if [[ "$1" == "NIH" ]] ; then
          shift ; notsource="EBI-Alt"

        elif [[ "$1" == "EBI" ]] ; then
          shift ; notsource="NIH-Alt"

        elif [[ "$1" == "local" ]] ; then
          shift ; localmode=true

        # Override default local_dir (dir source of the local files to be copied from)
        elif [ -d "$1" ]; then
          shift ; localmode=true ; local_dir="$1"

          if [[ "${local_dir}" == "${reflibdir}/genomes" ]] ; then
            echo "*** ERROR -- local_dir location cannot be the same as the reference/genomes location."
            valid=false
            break
          fi

        fi

        opt="opt5" ; continue
        ;;

      *)      # pseudo opt5 ; really any arg not recognized
        valid=false
        break
    esac
  done

  # Convert onefile to a valid genomes.csv file index
  if $valid && [[ -n $onefile ]] && ! $all ; then    # convert file name to numeric index
    found=false
    for ((j=1 ; j<genlen ; j++)); do
      if [[ "${onefile}" == "${finalf[$j]}" || "${onefile}" == "${initf[$j]}" ]] &&
         [[ $notsource != "${source[$j]}" ]] ; then
        oneindex=$j     # Save index entry for found file and break
        found=true
        break
      fi
    done

    if ! $found ; then
      echo "*** ERROR: $onefile is not a valid file name entry in the genomes.csv file."
      valid=false
    fi
  fi

  # If arguments invalid, then give usage
  if ! $valid ; then
    usage
    exit 1
  fi

}


#---------------------------------------------------------------------------------------------------------------
# "Library" mode. Will loop and perform actions queried from the user. Main one is to get and process
# refgenomes selected from genomes.csv. Also can list and bulk delete the refgenome library content.

library_mode() {
  local refgen rg
  local -a option
  local -i i

  refgen="${reflibdir}/genomes"
  case $OSTYPE in
    msys* | cygwin*) refgen=$("${oWGSEFIN}/cygwin64/bin/cygpath.exe" -w "${refgen}")
  esac
  # Dump banner for Library calls
  echo
  echo --------------------------------------------------------------------------------
  echo "WGS Extract Reference Genome Library Manager v${libver}"
  echo --------------------------------------------------------------------------------
  echo "WGSE v5.${currentVer} ${currentDate}; library at ${refgen}"
  echo "(*) are recommended. \"3x\" are uncompressed. See the Users Manual for more info"
  echo

  declare -i menu_cnt=0         # Redisplay the menu after so many downloads; so keep count

  # Setup menu and loop in BASH; majority of entries are from genomes.csv file menopt array

  printf -v PS3 '\nChoose which Reference Genome(s) to process now (1 to Exit): '
  option=( "Exit" "List Reference Genome Library" "Verify Reference Genome Library"
          "Recommended (*) (@US NIH)" "Recommended (*) (@EU EBI)"
           "${menopt[@]:1}" "Delete Reference Genome(s)")
  #  adjust=5   <--- change at top of file if changing # of "option"s before genomes list here

  select rg in "${option[@]}"; do
    case $rg in
      "Exit")
        break
        ;;

      "List Reference Genome Library")
        list_refgenomes
        menu_cnt+=1     # Do not want to scroll the list off with a new menu
        ;;

      "Verify Reference Genome Library")
        verify_refgenomes
        ;;

      "Recommended (*) (@US NIH)")
        for ((i=1 ; i<genlen ; i++)); do    # Go through all the menu options added by genomes.csv
          # These strings MUST match the menopt entries in genomes.csv EXACTLY (order does not matter)
          if [[ "hs38 (Nebula) (NIH) (*)" == "${menopt[$i]}" ||
                "hs37d5 (Dante) (NIH) (*)" == "${menopt[$i]}" ||
                "T2T_v2.0 (*)" == "${menopt[$i]}" ]]; then
            handle_refgenomes "$i"
          fi
        done
        menu_cnt+=3
        echo
        echo "Finished with Recommended (@US NIH)."
        ;;

      "Recommended (*) (@EU EBI)")
        for ((i=1 ; i<genlen ; i++)); do    # Go through all the menu options added by genomes.csv
          # These strings MUST match the menopt entries in genomes.csv EXACTLY (order does not matter)
          if [[ "hs38 (Nebula) (EBI) (*)" == "${menopt[$i]}" ||
                "hs37d5 (Dante) (EBI) (*)" == "${menopt[$i]}" ||
                "T2T_v2.0 (*)" == "${menopt[$i]}" ]]; then
            handle_refgenomes "$i"
          fi
        done
        menu_cnt+=3
        echo
        echo "Finished with Recommended (@EU EBI)."
        ;;

      "Delete Reference Genome(s)")
        # Do everything we can to make sure directory exists, is not the root (more than 2 characters), etc
        if [[ -d "${reflibdir}genomes/" ]] && (( ${#reflibdir} >= 3 )); then
          read -p "Delete which reference genome [Menu index, 0 for all]? " -n 3 -r # now requiring newline so no echo
          if [[ $REPLY =~ ^[0-9]+$ ]] ; then

            # Delete ALL reference genomes in the library
            if (( REPLY == 0 )) ; then
              ${bashx} "$process_refgenomes" "$debugstr" "clall"    # Asks for confirmation

            # Delete just the specified reference genome from the library given the menu index
            else
              index=$(( REPLY - adjust ))     # Convert Library menu ID to genomes.csv ID
              if is_index_good "$index" "library_delete" ; then

                # In case an unprocessed file, need to use the initf entry (if same as final then simply does final)
                if [[ -f "${initf[$index]}" ]] ; then
                  # We presume no finalf processed values exist to be deleted
                  ${bashx} "$process_refgenomes" "$debugstr" "clall" "${initf[$index]}"

                elif [[ -f "${finalf[$index]}" ]] ; then
                  ${bashx} "$process_refgenomes" "$debugstr" "clall" "${finalf[$index]}"

                else
                  echo "*** WARNING: nothing to delete for index $REPLY"
                fi

              # else  -- already reported error of a bad index in is_index_good
              fi
            fi

          else
            echo "*** ERROR: specified a non-numeric entry"
          fi
        else
          echo "*** ERROR: valid genomes library not found at ${reflibdir}genomes/"
        fi
        menu_cnt+=1
        ;;

      *)    # Main menu is not static (read from genomes.csv) so need to loop through list to find match (if at all)
        found=false
        for ((i=1 ; i<genlen ; i++)); do          # Go through all the menu options added by genomes.csv
          if [ "$rg" == "${menopt[$i]}" ]; then   # If we find the match, use the index to call the script
            handle_refgenomes "$i"                # Ignore return value as error already printed; library mode
            menu_cnt+=1
            found=true
            break
          fi
        done

        if ! "$found" ; then
          echo
          echo "Please enter a valid option (emter 1 to exit)"
        fi
        ;;
    esac

    echo

    if (( menu_cnt > 2 )); then
      # Cause the menu to be redisplayed again before the prompt
      REPLY=""
      menu_cnt=0
    fi

  done
  echo "Finished installing and processing Reference Genomes."
}


#------------------------------------------- ACTUAL MAIN CODE HERE  ------------------------------------------------

# Find Reflib and cd there; defined in zcommon.sh because needed everywhere
find_reflibdir  # Sets reflib with value from WGSE user settings or the default installation location if not set
# Must work out of reflibdir/genomes for clean and clean_all commands
\cd "${reflibdir}"genomes/ || { echo "*** ERROR: Internal, cannot cd to ${reflibdir}genomes/" ; exit ; }

# Get current release info for reflib and program (to use in banner)
read_current_release_info reflib "${reflibdir}reflib.json" > /dev/null   # Supress found message
libver=${currentVer}
read_current_release_info program "${oWGSEFIN}/program/program.json" > /dev/null   # Supress found message

# Read in genomes.csv file and transpose rows to columns; sets one array variable per column
read_genomes_file
genlen=${#initf[@]}   # Skipping first entry 0; so range of valid indices is 1 to (genlen-1)

# WGSE environment setup; so lets process the arguments. If it fails, error / usage is already printed and exits bad
process_args "$@"

# If debug mode on, echo the setup based on the arguments
if $debugmode ; then
  mesg="*** DEBUGMODE ON ***"
  [[ " ${stages[*]} " =~ " get " ]] && mesg+=" Get"
  [[ " ${stages[*]} " =~ " gver " ]] && { [[ -z $mesg ]] && mesg+=" Gver" || mesg+=", Gver" ; }
  [[ " ${stages[*]} " =~ " proc " ]] && { [[ -z $mesg ]] && mesg+=" Proc" || mesg+=", Proc" ; }
  [[ " ${stages[*]} " =~ " pver " ]] && { [[ -z $mesg ]] && mesg+=" Pver" || mesg+=", Pver" ; }
  (( ${#mesg} < 5 )) && mesg+=" Only"
  if ! $library ; then
    $stagecentric && mesg+=" *** Stage-Centric" || mesg+=" *** File-Centric"
    (( oneindex )) && mesg+=" *** File Index: ${oneindex}" || mesg+=" *** ALL Files"
  fi
  $localmode && mesg+=" *** LOCAL" || mesg+=" *** URL"
  echo "${mesg} ***"
fi

# Call either user interactive library mode or command line mode
if $library ; then
  library_mode                        # Loop through user GUI selections and make command line calls appropriately
else
  handle_refgenomes "$oneindex"       # If $all is true; $oneindex is 0 but does not matter as it is ignored
fi
exit $?
