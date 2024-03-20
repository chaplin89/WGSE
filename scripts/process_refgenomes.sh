#!/usr/bin/env bash
# Process Reference Genomes
# A standalone or sourced script to process Reference Genome files or folder for use in Bioinformatic tools
#  Generally called from get_refgenomes.sh overall script implementing the Library command, etc
# 
# If a directory is specified, the script naively presumes any file in there with the correct extension
#   is a FASTA file. It uses htsfile as a validity check after this first level screening. This is unlike
#   get_refgenomes which only works off the genomes.csv entries.
#
# This script does not use / depend on the genomes.csv file. Processes any file. Historically used current directory.
#    Was modified to allow path. But since restricted to only working in reflibdir.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
# Copyright (c) 2020-23 Randolph Harr
#

# -----------------------------------------------------------------------------------------------------------------
# Programs used (that must be available on the PATH):
# gunzip, unzip, bunzip2, 7z, sort, cut, (g)awk, sed, grep, tail, head, wc, md5sum, rm, bgzip, htsfile, samtools

np=16         # Number of Processors to use when available for bgzip  (Todo should read from system)
LANG=POSIX		# Needed for sort command

usage() {
  printf "Usage: %s [debug] {option}   where option is one of:\n" "$1"
  echo
  printf "      <directory>   # process ALL FASTA files in a directory; use dot (.) for current.\n"
  printf "   or <file(s)>     # List of 1 or more FASTA files to process.  Can have a path spec.\n"
  printf "   or clean [file]  # Remove all files created by this  script  in  the  ref  library.\n"
  printf "   or clall [file]  # Clean created AND reference genome files  in  the  ref  library.\n"
  echo
  printf " Files with extension fa, fna, fasta, gz, zip, bz, bz2, or 7z are checked  if  FASTAs.\n"
  echo
  printf " Processes Reference Genome FASTA files to assure they have  the  correct  compression\n"
  printf "  format. Then creates the DICT file and indices needed  by  tools.  Finally  catalogs\n"
  printf "  the file for WGSE use.  Processing a  whole  directory  causes  the  catalog  to  be\n"
  printf "  cleared first (mini-clean). Otherwise, the catalog is appended  to  with  any  newly\n"
  printf "  processed files. Only non-existing or older-than-the-FASTA files are (re)created.\n"
  echo
  printf " DICT and Index files are created alongside the source FASTA.  WGSE  catalog  info  is\n"
  printf "  created in the specified directory or current if specifying files.  Use files in the\n"
  printf "  current directory and without path specifications to assure both are together.\n"
  echo
  printf " \"clean\" and \"clall\" (clean all) with a file spec cleans only that  reference  genomes\n"
  printf "  files created by this script (and downloaded reference genome if \"clall\").  Default is\n"
  printf "  to operate on all reference genome files in the reference library.\n"
  echo
}

# Find the installation directory (normally part of zcommon.sh but when run standalone, cannot find that script)
WGSEDIR=$(/usr/bin/dirname "${BASH_SOURCE[0]}")   # Get the script location to determine install directory
WGSEABS=$(\cd "$WGSEDIR" || true; \pwd -P)                  # By cd'ing to it, resolve any aliases and symlinks
[[ $(/usr/bin/basename "${WGSEABS}") == "scripts" ]] && WGSEFIN=$(/usr/bin/dirname "${WGSEABS}") || WGSEFIN="${WGSEABS}"

declare pythonx
declare oWGSEFIN
declare reflibdir
declare rmfx
if ! source "${WGSEFIN}/scripts/zcommon.sh" dummy ; then
  echo "*** ERROR sourcing scripts/zcommon.sh in ${BASH_SOURCE[0]}. Exiting."
  exit 1
fi
cd "$reflibdir" || true   # Todo restore operating from anywhere and use reflibdir path if none specified and not .

#
# NOTE: We have commented out the BWA Index command due to the extensive resource usage (CPU and file space) for
#  doing BWA indices.  It is simply part of the main python code and run when needed.  Also different aligners have
#  different index files as it is.
#


# ------------------------------------- PROCESS COMMAND LINE PARAMETERS -------------------------------------------

# Process all parameters to see if expansion yields a list of one or more viable files to process.
# That way we know when to respond with usage. clean and clean_all work only in the current directory.
# The WGSE unique files are created in the current directory. Otherwise, DICT and index files are
# created with the reference genome source file. Relative / fixed paths are allowed for dir and files.

shopt -s nullglob
sumFB="WGSE"        # Main file base to be built on for WGSE catalog files
sumFPB="./WGSE"     # Full path with base name to build on (when a directory is specified; default is current dir)
debugmode=false

if [[ "$1" == "debug" ]] ; then
  script=$(basename "$0")
  echo "***DEBUG: debug mode turned on in $script"
  shift
  debugmode=true

elif [[ -z "$1" ]] ; then   # Happens when called with "" in place of debug
  shift
fi

# Clean files created by this script and return
if [[ "$1" == "clean" ]] ; then
  $debugmode && echo "***Performing Clean $2"

  # Single reference genome specified
  if [ -f "$2" ] ; then
    base=${2%%.*}   # Get base name for special files
    $rmfx "$2".fai "$2".gzi "$base".dict "$base"_n*.csv "$base".wgse ./HTS.tmp || true   # This scripts files
    $rmfx "$2".amb "$2".ann "$2".bwt "$2".pac "$2".sa || true                            # BWA Index files

  # Whole directory (all related reference genome files created by this script; left only with downloaded ref genomes)
  else                    # Whole directory operation
    $rmfx ./*fai ./*gzi ./*dict ./*_nbin.csv ./*_ncnt.csv ./HTS.tmp || true              # This scripts files
    $rmfx ./*amb ./*ann ./*bwt ./*pac ./*sa || true                                      # BWA Index files
    $rmfx ./*wgse "${sumFPB}.csv" "${sumFPB}"_dict*.csv "${sumFPB}"_uniq*.csv || true    # WGSE library files

  fi

  (return 0 2>/dev/null) && return 0 || exit 0

# Clean files created here, BWA indices if here also, AND the reference genomes themselves; and return
elif [[ "$1" == "clall" ]] ; then
  $debugmode && echo "***Performing Clean_all (clall) $2"

  # Single reference genome operation
  if [[ -f "$2" ]] ; then
    base=${2%%.*}   # Get base name for special files
    $rmfx "$2".fai "$2".gzi "$base".dict "$base"_n*.csv "$base".wgse ./HTS.tmp || true   # This scripts files
    $rmfx "$2".amb "$2".ann "$2".bwt "$2".pac "$2".sa || true                            # BWA Index files
    $rmfx "$2" || true      # todo pre-process and post-process file name                # Actual Reference Genome

  # Whole directory (all related reference genomes files AND downloaded reference genomes)
  else
    $rmfx ./*fai ./*gzi ./*dict ./*_nbin.csv ./*_ncnt.csv ./HTS.tmp || true              # This scripts files
    $rmfx ./*amb ./*ann ./*bwt ./*pac ./*sa || true                                      # BWA Index files
    $rmfx ./*wgse "${sumFPB}.csv" "${sumFPB}"_dict*.csv "${sumFPB}"_uniq*.csv || true    # WGSE library files

    read -p "Do you really want to remove ALL the reference genomes in $PWD [y/N]? " -n 1 -r ; echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      echo ' ... deleted all downloaded reference genome files in the reference library.'
      $rmfx ./*fa.gz ./*fasta.gz ./*fna.gz || true    # Delete processed reference genomes
      $rmfx ./*.fa ./*.fasta ./*.fna ./*.zip || true  # in case any unprocessed (initfn) also

    else
      echo ' ... leaving reference genomes alone.'
    fi
  fi

  (return 0 2>/dev/null) && return 0 || exit 0

# else process the single directory specified ; most commonly will be "."
elif [ -d "$1" ] ; then                   # Create a file list from directory content
  file_list=("${1%/}"/*.{fa,fna,fasta,7z,gz,zip,bz,bz2})    # Pattern match for known extensions of reference genomes.

  if [[ "$1" != "." ]] ; then
    sumFPB="${1%/}"/${sumFB}      # Removes any possible trailing slash before adding one back in
  fi

  # Want to process whole directory; so start with a clean slate of summary files (catalog cleared)
  # Different than earlier because now we have to look into a specific directory; not necessarily the current (".")
  $rmfx "${sumFPB}.csv" "${sumFPB}_dict.csv" "$1"/*wgse || true
  $rmfx "${sumFPB}"_uniq*.csv "${sumFPB}"_dict*.csv || true

  wholedir=true

# else must be a list of one or more files; so create a list of files from the parameters
elif (( $# > 0 )) ; then
  file_list=("$@")
  wholedir=false

fi

# Now that we have processed the parameters, we can determine usage (in case no files listed after processing)
if (( $# == 0 || ${#file_list[@]} == 0 )) ; then
  script=$(basename "$0")
  usage "$script"
  (return 0 2>/dev/null) && return 1 || exit 1
fi


# ----------------------------------- PRIMARY CHROMOSOME AND MODEL DEFINITIONS ----------------------------------
# Setup variables for finding primary chromosomes.  Use LN entries as main source. But some alt contigs have
#  same length.  So also then filter by common SN forms. Should change to read in file / table more easily updated
#
# Build   15/(33)        16/34         17/35         18/36         19/37        (20)/38     hg01243 v3    chm13 v0.9    chm13 v1.1    chm13 v1.0 (chm13 v0.7 X, hg002 v2 XY, hg002 v2.7 XY)
chr[0]="LN:245203898\|LN:246127941\|LN:245522847\|LN:247249719\|LN:249250621\|LN:248956422\|LN:248415701\|LN:248387561\|LN:248387328\|LN:248387497"  # chr1
chr[1]="LN:243315028\|LN:243615958\|LN:243018229\|LN:242951149\|LN:243199373\|LN:242193529\|LN:242509959\|LN:242696759\|LN:242696752\|LN:242696747"  # chr2
chr[2]="LN:199411731\|LN:199344050\|LN:199505740\|LN:199501827\|LN:198022430\|LN:198295559\|LN:200717518\|LN:201106621\|LN:201105948\|LN:201106605"  # chr3
chr[3]="LN:191610523\|LN:191731959\|LN:191411218\|LN:191273063\|LN:191154276\|LN:190214555\|LN:193408891\|LN:193575384\|LN:193574945\|LN:193575430"  # chr4
chr[4]="LN:180967295\|LN:181034922\|LN:180857866\|LN:180857866\|LN:180915260\|LN:181538259\|LN:182049998\|LN:182045443\|LN:182045439\|LN:182045437"  # chr5
chr[5]="LN:170740541\|LN:170914576\|LN:170975699\|LN:170899992\|LN:171115067\|LN:170805979\|LN:171893897\|LN:172126875\|LN:172126628\|LN:172126870"  # chr6
chr[6]="LN:158431299\|LN:158545518\|LN:158628139\|LN:158821424\|LN:159138663\|LN:159345973\|LN:160394084\|LN:160567465\|LN:160567428\|LN:160567423"  # chr7
chr[7]="LN:145908738\|LN:146308819\|LN:146274826\|LN:146274826\|LN:146364022\|LN:145138636\|LN:146097661\|LN:146259347\|LN:146259331\|LN:146259322"  # chr8
chr[8]="LN:134505819\|LN:136372045\|LN:138429268\|LN:140273252\|LN:141213431\|LN:138394717\|LN:149697505\|LN:150617238\|LN:150617247\|LN:150617274"  # chr9
chr[9]="LN:135480874\|LN:135037215\|LN:135413628\|LN:135374737\|LN:135534747\|LN:133797422\|LN:134341430\|LN:134758139\|LN:134758134\|LN:134758122"  # chr10
chr[10]="LN:134978784\|LN:134482954\|LN:134452384\|LN:134452384\|LN:135006516\|LN:135086622\|LN:134654341\|LN:135129789\|LN:135127769\|LN:135127772" # chr11
chr[11]="LN:133464434\|LN:132078379\|LN:132449811\|LN:132349534\|LN:133851895\|LN:133275309\|LN:133439878\|LN:133324792\|LN:133324548\|LN:133324781" # chr12
chr[12]="LN:114151656\|LN:113042980\|LN:114142980\|LN:114142980\|LN:115169878\|LN:114364328\|LN:113815969\|LN:114240132\|LN:113566686\|LN:114240146" # chr13
chr[13]="LN:105311216\|LN:105311216\|LN:106368585\|LN:106368585\|LN:107349540\|LN:107043718\|LN:100860689\|LN:101219190\|LN:101161492\|LN:101219177" # chr14
chr[14]="LN:100114055\|LN:100256656\|LN:100338915\|LN:100338915\|LN:102531392\|LN:101991189\|LN:99808683\|LN:100338336\|LN:99753195\|LN:100338308"   # chr15
chr[15]="LN:89995999\|LN:90041932\|LN:88827254\|LN:88827254\|LN:90354753\|LN:90338345\|LN:96296229\|LN:96330509\|LN:96330374\|LN:96330493"           # chr16
chr[16]="LN:81691216\|LN:81860266\|LN:78774742\|LN:78774742\|LN:81195210\|LN:83257441\|LN:83946371\|LN:84277212\|LN:84276897\|LN:84277185"           # chr17
chr[17]="LN:77753510\|LN:76115139\|LN:76117153\|LN:76117153\|LN:78077248\|LN:80373285\|LN:80696073\|LN:80537682\|LN:80542538\|LN:80542536"           # chr18
chr[18]="LN:63790860\|LN:63811651\|LN:63811651\|LN:63811651\|LN:59128983\|LN:58617616\|LN:61612450\|LN:61707413\|LN:61707364\|LN:61707359"           # chr19
chr[19]="LN:63644868\|LN:63741868\|LN:62435964\|LN:62435964\|LN:63025520\|LN:64444167\|LN:67262993\|LN:66210261\|LN:66210255\|LN:66210247"           # chr20
chr[20]="LN:46976537\|LN:46976097\|LN:46944323\|LN:46944323\|LN:48129895\|LN:46709983\|LN:44996062\|LN:45827694\|LN:45090682\|LN:45827691"           # chr21
chr[21]="LN:49476972\|LN:49396972\|LN:49554710\|LN:49691432\|LN:51304566\|LN:50818468\|LN:51228122\|LN:51353916\|LN:51324926\|LN:51353906"           # chr22
chr[22]="LN:152634166\|LN:153692391\|LN:154824264\|LN:154913754\|LN:155270560\|LN:156040895\|LN:154343774\|LN:154259664\|LN:154259566\|LN:154259625\|LN:154269076\|LN:154349815\|LN:154434329"  # chrX
chr[23]="LN:50961097\|LN:50286555\|LN:57701691\|LN:57772954\|LN:59373566\|LN:57227415\|LN:62480187\|LN:62456832\|LN:62460029"                        # chrY (no CHM13 entries)
# Need separate entries later so define array first before concatinating to form chrln match string here
chrln=${chr[0]}
for val in "${chr[@]:1}" ; do
  chrln+="\|${val}"
done

# Some alt contigs have same length as a chromosome; so further refine selection by known primary entry SNs (prefix)
chrsnhg="SN:CHR[1-9XY]\>\|SN:CHR1[0-9]\>\|SN:CHR2[0-2]\>\|SN:CHRX_V0.7\>"  # HGP naming (special CHM13 v0.7)
chrsneb="SN:[1-9XY]\>\|SN:1[0-9]\>\|SN:2[012]\>"                           # EBI naming
chrsnnc="SN:NC_0000[0-2][0-9]"                                             # RefSeq naming (Accession)
chrsngb="CM0006\|CM0349[567]\|CP0682[567]\|CP08656[89]\|SN:CHR[XY]_HG002"  # Genbank naming (Accession)
chrsn="${chrsnhg}\|${chrsneb}\|${chrsnnc}\|${chrsngb}"
# Alt contig names in HGP and EBI can start with full chromosome name; so make sure checking whole field (>)
# In accession naming, alt contigs have their own accession entries that start differently than the chromosomes

# The above is only for the primary (chromosomes).  Need to process the mitochondrial model often there also.
mitln="LN:16569\>\|LN:16571\>\|LN:16568\>"
mitsn="SN:CHRM\>\|SN:CHRMT\>\|SN:MT\>\|[:|]J01415\|SN:CP068254\|SN:CM032116\|SN:NC_001807\|SN:NC_012920"


# --------------------------------- MAIN PROCESSING START ---------------------------------------------------------
# Start processing ... first, setup global files if they do not exist.  Then loop on each argument / refgen file

# Now that we are starting; start tracking errors in file processing
declare -i exit_status=0

if $wholedir ; then
  printf "================================================================================\n"
  printf "Processing Reference Genome file(s) ...\n"
  $debugmode && echo "***DEBUG: ${#file_list[@]} parameter(s): ${file_list[*]}"
else
  printf "Processing Reference Genome file %s \n" "${file_list[@]}"
fi

# Setup initial TSV summary file header if not yet written; for when initial processing of whole directory
# Columns are: File Name, Major/Minor build code, SN Cnt ;
#              MD5Sums: SN/LN, SN/LN/M5, LN/M5, Chromo LN/M5, Chromo SN, initial File, final file ;
#              Mito: SN, LN, M5 ; Error Message(s) ; list of 24 primary SNs
# Chromo LN/M5 identifies Major Build, Full LN/M5 the Minor Class, SN Cnt the patch level
if [[ ! -f "${sumFPB}.csv" ]] ; then
  printf $'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "File" "Build" "SN_CNT" "BAM_(SN,_LN)" "CRAM_(SN,_LN,_M5)" "FASTA_(LN,_M5)" "Chromo_(LN,_M5)" "Chromo_(SN)" \
    "Orig_File_MD5" "Final_File_MD5" "Mito:_SN" "LN" "M5" "Message" "SNs" > "${sumFPB}.csv"
fi


# ------------------------------ LOOP THROUGH EACH FILE TO PROCESS --------------------------------------------------
#
$debugmode && echo ""
for rgfile in "${file_list[@]}" ; do
  $debugmode && echo "***DEBUG: file is $rgfile ; cwd is $PWD"

  # For common case of using current directory, strip off ./ that starts each file name
  if [[ "${rgfile:0:2}" == "./" ]] ; then
    rgfile=${rgfile:2}
  fi

  # fpa=${rgfile%/*}    # File path (if it exists)  (not used; but keep in case we need it later)
  fbn=${rgfile##*/}   # File base name
  ext=${rgfile##*.}   # File (last) extension

  # Check if even a file (in case am imprecise glob match / directory spec pulled in non-files to just skip)
  if [[ ! -f "$rgfile" ]] ; then        # Always skip non-files
    $debugmode && echo "${fbn}: Skipping file that does not exist or is not a file"
    if [[ ! $debugmode && ! $wholedir ]] ; then         # Report an error IFF gave explicit file (list); not directory
      echo "${fbn}: ***ERROR: Skipping a non-file"
      exit_status+=1

    fi
    continue

  fi

  # Check if file content and extension correct
  # here-documents from variable (<<< $HTS) not working in subsetted cygwin bootstrap; so change to temp file
  # Instead of HTS.tmp file, use HTS=$(htsfile "$file") and <<< $HTS if here-documents working
  htsfile "$rgfile" > "HTS.tmp"

  if [[ "$ext" =~ ^(zip|7z|bz|bz2)$ ]]; then
    # htsfile cannot see in (to tell if FASTA) for .zip, .7z nor .bz2 archive formats (but does recognize bzip2 format)
    # So go ahead if those and check if a FASTA after recompressing; might recompress some non-bioinformatic files but ...
    # that is a risk the user takes if they do a general directory specification
    echo "${fbn}: ***WARNING: Cannot look into .zip/.bz/.bz2/.7z files; so recompressing to peak inside."

  elif grep -v "FASTA" < "HTS.tmp" &> /dev/null ; then
    echo "${fbn}: ***ERROR: Skipping a non-FASTA file"
    $rmfx "HTS.tmp" || true
    exit_status+=1
    continue

  fi  # simply fall through as next check determines if you need to recompress or not

  md5i=$(md5sum "$rgfile" | cut -d' ' -f1)    # Save initial file md5sum

  # ---------------------------------------------------------------------------------------------------------------
  # Check for bgzip compression form; (re)compress if not in that form
  if grep "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
    $debugmode && echo "   $rgfile: Already a BGZF compressed FASTA file!"   # simply keep silent
    filen="$rgfile"

  else
    # OK, if here, is (likely) a FASTA but is not BGZF compressed
    case "$ext" in
      gz)     form="gzip"    ;;
      zip)    form="pkzip"   ;;
      7z)     form="7zip"    ;;
      bz)     form="bzip"    ;;
      bz2)    form="bzip2"   ;;
      fa)     form="none"    ;;
      fna)    form="none"    ;;
      fasta)  form="none"    ;;
      *)		# Well that is embarrassing
	      echo "$fbn: ***ERROR: Unknown file extension ($ext) (internal error)"
	      $rmfx "HTS.tmp" || true
	      exit_status+=1
	      continue
	      ;;
    esac
    echo "   $fbn: Fixing Reference Genome File compression ($form to BGZF)"
    
    case "$ext" in
      gz)		# BGZF already checked above; so must be gzipped
        { gzip -d -f "$rgfile" ; bgzip -if@ $np "${rgfile%.*}" ; } && filen="$rgfile" ;;
  	    # This is the only format where you end up with the same file extension; they should have used .bgz but ...

      zip)		# Likely zip'ped 
	      { unzip -p "$rgfile" | bgzip -cf@ $np > "${rgfile%.*}.gz" ; } && $rmfx "$rgfile" && filen="${rgfile%.*}.gz"  ;;

      7z)		# A popular, very high compression format not normally seen in Bioinformatics
	      { 7z e -mmt$np -so "$rgfile" | bgzip -cf@ $np > "${rgfile%.*}.gz" ; } && $rmfx "$rgfile" && filen="${rgfile%.*}.gz"  ;;

      bz | bz2)		# A popular Linux/Unix/MacOSX high compression format sometimes seen in Bioinformatics
        { bunzip2 -cf "$rgfile" | bgzip -cf@ $np > "${rgfile%.*}.gz" ;} && $rmfx "$rgfile" && filen="${rgfile%.*}.gz"  ;;

      fa | fna | fasta)	# Appears likely not compressed; but will check
  	    if grep -- "-compressed" < "HTS.tmp" &> /dev/null ; then
	        echo "$fbn: ***ERROR: File is compressed but with only a FASTA extension (no .gz extension). Skipping. "
	        $rmfx "HTS.tmp" || true
	        exit_status+=1
	        continue
	      fi
  	    bgzip -if@ $np "$rgfile"
	      filen="$rgfile.gz"
  	    ;;
    esac
    # Above does not handle .zip / .7z archive where multiple files (i.e. like tar+gzip combo from original pkzip)
    #  Singular, internal file name must be the same name as archive file without extension

    # Check again if matched file with extension is a valid, compressed FASTA file just to make sure
    if [[ ! -f "$filen" ]] ; then
      echo "$fbn: ***ERROR: Compression conversion failed (missing file). Skipping"
      $rmfx "HTS.tmp" || true
      exit_status+=1
      continue
    fi
    #HTS=$(htsfile "$filen")
    htsfile "$filen" > "HTS.tmp"
    if grep -v "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
      echo "$fbn: ***ERROR: Compression conversion failed (not BGZF). Skipping"
      $rmfx "HTS.tmp" || true
      exit_status+=1
      continue
    fi
  fi

  md5fc=$(md5sum "$filen" | cut -d' ' -f1)
  $rmfx "HTS.tmp" || true

  # ---------------------------------------------------------------------------------------------------------------
  # Check and (re)create any supporting indices
  #
  # By this point, Ref Genome file ends in .gz, is a FASTA and is BGZF compressed. Does allow explicitely named file
  #  in parameter list that is a FASTA and BGZF compressed but without the proper extension (e.g. .fa.gz).
  #  Some programs will not allow that but who are we to judge.
  #
  # If file was properly compressed AND files to be created are up to date, nothing will get posted about this file!

  # (Re)create any missing or old indices (using a new filename $filen if compression changed)
  # GATK uses the .dict file and wants it to be named without {fasta,fa,fna}.gz extension so we oblige
  # Reminder: -nt is the BASH "file newer than" operator . File name may have more dots than just our extensions.
  filed=$(echo "$filen" | sed "s/.fasta.gz//;s/.fna.gz//;s/.fa.gz//")  # strip known extension to add .dict back on
  fileb=${filen##*/}         # We use the whole filename for status reporting; so if changed need a new basename
  if [ "$filen" -nt "${filed}.dict" ]; then
    echo "   $fileb: Creating FA DICTionary"
    samtools dict "$filen" -o "$filed.dict"
    (( $? != 0 )) && { exit_status+=1 ; continue ; }
  fi
  if [ "$filen" -nt "${filen}.fai"  ]; then
    echo "   $fileb: Creating FA Index (FAI)"
    samtools faidx "$filen"    # Usually creates .gzi as well
    exit_status+=$?
  fi
  if [ "$filen" -nt "${filen}.gzi"  ]; then
    echo "   $fileb: Creating BGZip Index (GZI)"
    bgzip -r "$filen"          # Also samtools index works
    (( $? != 0 )) && { exit_status+=1 ; continue ; }
  fi
  if [ "$filen" -nt "${filed}_ncnt.csv" ]; then
    echo "   $fileb: Counting N's (ncnt, nbin, bed)"
    case $OSTYPE in
      msys* | cygwin*)   ofilen=$(cygpath -m "${filen}")  ;;
      darwin* | linux*)  ofilen="${filen}"                ;;
    esac
    "${pythonx}" "${oWGSEFIN}/program/countingNs.py" "${ofilen}" | grep -v "INFO"  # Grep returns 1 if nothing to filter
    (( PIPESTATUS[0] != 0 )) && { exit_status+=1 ; continue ; }
  fi
  # BWA Index creates .bwt (30 min, 3 GB), .pac (800MB), .ann, .amb, and .sa (10 min, 1.5GB); 6.4GB total!
  # [ "$filen" -nt "$filen.bwt" ]  && echo "$filen: Creating BWA Indices: 3 hours, 5.5 GB" && bwa index "$filen"
  # Todo add parameter to optionally turn on BWA index capability

  # ------------------------------ CREATING WGS EXTRACT CATALOG INFO ------------------------------------------------
  # from the DICT file; when processing a directory we cleared all the .wgse files first
  [ "$filen" -nt "$filed.wgse" ] && echo "   $fileb: Creating WGS Extract Catalog Info"  && {

    # ISSUE: How to best sort DICT file to get matching md5 hash on primary chromosomes of the same minor class?
    #   We numeric sort on LN only now. So md5 hash for same minor class primary chromosomes matches. And can recreate
    #   LN md5 with the BAM. We started with a simple "sort -d" on whole line (as preferred -V was not available
    #   on all platforms) but that uses the SN field as the primary sort order. So md5 hashes of same minor class
    #   but different SN naming appear different. Sorting by LN descending gets closer to what visually you would
    #   like to see of the primary chromosomes roughly in order.  But it still allows some contigs (patches) that
    #   are the same length as the primary chromosome to be interspersed.  And pushes the Mito sequence near the end.
    #
    # sort -d is how we did it in the begining. But that is an alphabetic sort. (e.g. chr2 comes after Chr10-19)
    # sort -V may be more optimal but -V is not portable (GNU Util only); and brings chrN_xxxx contigs in the middle
    # sort --key=2.4brn,2 -k1.4b,1 -k3.4b,3   # to sort by descending LN, then SN then M5 so minor class have same md5
    # Often no sort (original order) has all the primary chromosomes and mito first before contigs. But not all
    #   files of the same minor class have the same sort order.

    # <<< $var here-document using variable is not working in subsetted cygwin64 of Windows; so we use tmp file instead

    # ---------------------------------------------------------------------------------------------------------------
    # Calculate various md5 hashes from various forms of DICT subset

    # Throw out DICT first line, take only columns 2 to 4. Store temporarily. This is the key file.
    tail -n +2 "$filed.dict" | cut -f2-4 > "$filed.otmp"

    # Original order and case as in file. md5 hash changes if you sort and/or upcase first
    md5so=$(cut -f1     "$filed.otmp" | md5sum | cut -d' ' -f1)     # MD5Sum SN only (SN        )
    md5lo=$(cut -f2     "$filed.otmp" | md5sum | cut -d' ' -f1)     # MD5Sum LN Only (    LN    )
    md5bo=$(cut -f1-2 < "$filed.otmp" | md5sum | cut -d' ' -f1)     # MD5Sum BAM     (SN, LN    ) => BAM determinable
    md5co=$(cut -f1-3 < "$filed.otmp" | md5sum | cut -d' ' -f1)     # MD5Sum CRAM    (SN, LN, M5)
    md5fo=$(cut -f2-3 < "$filed.otmp" | md5sum | cut -d' ' -f1)     # MD5Sum FASTA   (    LN, M5)

    # Now upcase and sort to remove dependency on that. (Always the same?)
    awk '{print toupper($0)}' < "$filed.otmp" | sort -d > "$filed.stmp"

    # These are the original values calculated (upcase and sorted) from the first year (sort -d though)
    md5ss=$(cut -f1     "$filed.stmp" | md5sum | cut -d' ' -f1)     # MD5Sum SN only (SN        )
    md5ls=$(cut -f2     "$filed.stmp" | md5sum | cut -d' ' -f1)     # MD5Sum LN Only (    LN    )
    md5bs=$(cut -f1-2 < "$filed.stmp" | md5sum | cut -d' ' -f1)     # MD5Sum BAM     (SN, LN    ) => BAM determinable
    md5cs=$(cut -f1-3 < "$filed.stmp" | md5sum | cut -d' ' -f1)     # MD5Sum CRAM    (SN, LN, M5)
    md5fs=$(cut -f2-3 < "$filed.stmp" | md5sum | cut -d' ' -f1)     # MD5Sum FASTA   (    LN, M5)

    # Capture only Primary (Chromosomes) DICT entries ; filter by LN sizes then SN names (Orig not upcase so cannot use)
    grep -w -e "$chrln" "$filed.stmp" | grep -e "$chrsn"  > "$filed.pstmp"  # sorted order
    # Note: above retains multiple chrY values if they exist

    # Primaries only; sorted, upcase.  Common naming = md5sp. Common major build = md5lp, Common minor class = md5fp.
    md5sp=$(cut -f1   "$filed.pstmp" | md5sum | cut -d' ' -f1)  # MD5Sum Primary SN only (SN        ) (sort, up)
    md5lp=$(cut -f2   "$filed.pstmp" | md5sum | cut -d' ' -f1)  # MD5Sum Primary LN only (    LN    ) (sort, up)
    md5fp=$(cut -f2-3 "$filed.pstmp" | md5sum | cut -d' ' -f1)  # MD5Sum Primary FASTA   (    LN, M5) (sort, up)

    # "fixed" primaries: pull each chromosome in fixed order; sorted and upcase and take only first occurence
    for val in "${chr[@]}" ; do
      sn=$(grep -w -e "${val}" "$filed.stmp" | grep -e "$chrsn" | head -n 1)
      printf "%s\n" "${sn}" >> "${filed}.fptmp"
    done

    # Fixed order, primaries only.  Common naming = md5sf. Common major build = md5lp, Common minor = md5ff.
    md5sf=$(cut -f1   "$filed.fptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only SN only (SN        ) (fixed order)
    md5lf=$(cut -f2   "$filed.fptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only LN only (    LN    ) (fixed order)
    md5ff=$(cut -f2-3 "$filed.fptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only FASTA   (    LN, M5) (fixed order)

    # Mitochondrial entries (not captured in primary earlier); use sorted and upcased as multiple and take first later
    # Models are rCRS, Yoruba, and individual T2T / HPP entries for samples (CHM13, hg01243, etc)
    grep -w -e "$mitln" "$filed.stmp" | grep -e "$mitsn" > "$filed.mtmp"	# Keep all columns for detail

    # ---------------------------------------------------------------------------------------------------------------
    # Determine stats for file(s) and report on any errors

    # SN Count (whole DICT), Primary SN count (just chromosomes and mito), Y chromosome count, Mitochondria count
    declare -i snct pcnt ycnt mcnt
    errp=""                                             # Initialze error string so can always concatenate later

    # Calculate counts of various sequences ; redirect stdin on wc so filename not included
    snct=$(wc -l < "$filed.otmp")                       # True SN count; on full file of SN entries
    pcnt=$(wc -l < "$filed.pstmp")                      # Primary count; redirect stdin so file name not included
    ycnt=$(grep -c -w -e "${chr[23]}" "$filed.pstmp")   # Y chromosome count (sometimes 0, sometimes 2)
    mcnt=$(wc -l < "$filed.mtmp")   # Redirect stdin so filename not included.

    $debugmode && echo "***DEBUG: Primary Chromosomes: $pcnt" && head -1 "$filed.pstmp"

    # Determine what errors exist; save and/or report
	  if (( snct < 25 )); then    # Smallest Analysis model is 84 SN entries (24/25 in T2T)
	    echo "$fileb: ***WARNING: Too few SN entries (<25)***"
	    errp+=" ***WARN: Too few SN entries (<25)***"
	  fi

    # Only report one of these three errors as multiple may be common
    if (( pcnt == 0 )) ; then
      echo "$fileb: ***WARNING: No Chromosomes found; unrecognized LN lengths? Build previously seen?"
      errp+=" ***WARN: No Chromosomes found ***"

    elif (( ycnt != 1 )) ; then # First check if missing Y or more than one Y; report that as issue before checking pcnt
      # M and Y are only duplicates found so far; so give a special error if found to be Y here; M is later
      echo "$fileb: ***WARNING: 1 expected, $ycnt found: Y chromosome entries in ref model"
      errp+=" ***WARN:Y $ycnt!=1 ***"

    elif (( pcnt != 24 )) ; then
      echo "$fileb: ***WARNING: 24 expected, $pcnt found: chromosomes in primary ref model"
      errp+=" ***WARN:P $pcnt!=24 ***"

    fi

    if (( mcnt != 1 )) ; then
      echo "$fileb: ***WARNING: 1 expected, $mcnt found: mitrochondrial entries in ref model"
      errp+=" ***WARN:M $mcnt!=1 ***"
    fi
    if (( mcnt > 1 )) ; then     # Truncate to first entry to simplify looking later
      head -n 1 "$filed.mtmp" > "temp.mtmp" && mv "temp.mtmp" "$filed.mtmp"
    fi

    # ---------------------------------------------------------------------------------------------------------------
    # Determine Mito SN, LN and M5; use to determine SN type, mito build type (19 vs 37), etc.

    # Pull out just SN of Mito
    msn=$(cut -f1 "$filed.mtmp")
    chrMSN="${msn:3}"

    # Setup trailing ID based on Sequence Naming convention using Mito entry as primary identifier
    #  in case Mito does not exist, check for primary chromosome name entries using previous patterns
    case $chrMSN in
 	    CHRM)		     msn="h"   ;;  # UCSC chrN numeric / single-alphabetic naming (standard)
      CHRMT)		   msn="ht"      # UCSC chrN naming using oddball chrMT
        echo "$fileb: ***WARNING: chrMT name is non-standard ***"
        errp+=" ***WARN: chrMT name is non-standard ***"  ;;
      MT)          msn="g"   ;;  # EBI mumeric / single-alphabetic naming except MT
 	    NC_012920*)  msn="n"   ;;  # NCBI RefSeq Acquisition ID (rCRS)
 	    NC_001807*)  msn="n"   ;;  # NCBI RefSeq Acquisition ID (Yoruba)
 	    J01415.2)    msn="c"   ;;  # NCBI GenBank ID naming (T2T v2)
 	    CP068254.1)  msn="c"   ;;  # NCBI GenBank ID naming (T2T v1.x)
 	    CM032116.2)  msn="c"   ;;  # NCBI GenBank ID naming (T2T hg01243)
 	    "GI|113200490|GB|J01415.2|HUMMTCG")
                   msn="c"   ;;  # NCBI GenBank ID naming (T2T)
      *)  # Check if a ref model we already know about; maybe just no MT. Otherwise, truly unknown.
        if (( pcnt > 0 )) ; then    # Model we know about; so has to be one of the three types we already checked for
          declare -i SNhg SNeb SNnc
          SNhg=$(grep -c -e "$chrsnhg" "$filed.stmp")
          SNeb=$(grep -c -e "$chrsneb" "$filed.stmp")
          SNnc=$(grep -c -e "$chrsnnc" "$filed.stmp")
          if   (( SNhg > 0 )) ; then  msn="h"
          elif (( SNeb > 0 )) ; then  msn="g"
          elif (( SNnc > 0 )) ; then  msn="n"
          else                        msn="c"
          fi
        else
          msn="x"       # Unrecognized naming convention; or no MT in FASTA model file?
          unkSN=$(grep -e "SN\:" "$filed.stmp" | head -n 1)
          echo "$fileb: ***WARNING: Unregonized Sequence Naming (${unkSN})"
          errp+=" ***WARN: Unrecognized SN naming convention ***"
        fi  ;;
 	  esac

    chrMLN=$(cut -f2 "$filed.mtmp")     # Used in WGSE line print out later

   	# Determine mitochondria model Build variation based on mitochondrial model found
    chrMM5=$(cut -f3 "$filed.mtmp")
    case $chrMM5 in
      M5:C68F52674C9FB33AEF52DCF399755519)	mbuild="37"  ;; # rCRS build 37, 38, etc
      M5:D2ED829B8A1628D16CBEEE88E88E39EB)	mbuild="19"  ;; # Yoruba build 19, 18, etc
      M5:EC493A132AC4823AA696E37109F64972)  mbuild="99"  ;; # T2T model: chm13 v1.0, v1.1, v2.0 (and those that include it)
      M5:2AEA08C58600A30435E4302A82481DC0)  mbuild="99"  ;; # T2T model: chm13 v0.9
 	    M5:30F23EB261CAB50B34C0BED87EE38C7E)  mbuild="99"  ;; # T2T model: hg01243 v3
      *)                            	      mbuild="xx"     # Unrecognized model; RSRS?
        if [ -n "$chrMM5" ] ; then    # Only report error if set; otherwise no MT in model already reported
          echo "$fileb: ***WARNING: Unregonized Mitochondrial Model (${chrMM5})"
          errp+=" ***WARN: Unrecognized Mito Model***"
        fi ;;
    esac

    # ---------------------------------------------------------------------------------------------------------------
    # Determine major/minor Build type (ignoring SN names; based on md5 of primary chromosome LN and M5 fields only)
 	  case $md5fp in
      13cbd449292df5bd282ff5a21d7d0b8f)   build="T2Tv20a"       ;;  # T2T CHM13+HG002Y v2.0 (accession; diff sort order)
      1e34cdea361327b59b5e46aefd9c0a5e)   build="HG16"          ;;  # hg 16 / NCBI 34
      3566ee58361e920af956992d7f0124e6)   build="HG15"          ;;  # hg 15 / NCBI 33
      4136c29467b6757938849609bedd3996)   build="NCB38"         ;;  # NCBI 38 all patches (GENBANK, REFSEQ)
 	    46cf0768c13ec7862c065e45f58155bf)   build="EBI18"         ;;  # EBI 18
      4bdbf8a3761d0cd03b53a398b6da026d)	  build="HG38"  	      ;;  # HG38 all patches
      4bf6c704e4f8dd0d31a9bf305df63ed3)   build="THGv27"        ;;  # T2T CHM13 v1.1 with HG002 xy v2.7
      4d0aa9b8472b69f175d279a9ba8778a1)   build="HPPv11"        ;;  # HPP CHM13 v1.1 with GRCh38 Y
   	  591bb02c89ed438566ca68b077fee367)   build="1K37p"         ;;  # EBI GRCh37 Errored  (extra Y) (few are EBI37)
      5a23f5a85bd78221010561466907bf7d)   build="EBI37"         ;;  # EBI 37 (including hg 37 ENA / WGSE)
      5e16e3cbdcc7b69d21420c332deecd3b)   build="T2Tv10"        ;;  # T2T CHM13 v1.0 (original)
      5f451c1014248af62b41c18fec1c3660)   build="T2Tv07"        ;;  # T2T CHM13 v0.7 (chr X only; not a true build)
      65a05319ad475cf51c929d3b55341bc2)   build="THGv20"        ;;  # T2T CHM13 v1.1 with HG002 xy v2
      7083d4ee8aa126726961ab1ae41c66c1)   build="THG1243v3"     ;;  # T2T HG01243 (Pr1) (accession)
   	  7a5eb72fb45c4567431651aa6f9edfef) 	build="1K${mbuild}"   ;;  # 1K 19 / 37
      7cee777f1939f4028926017158ed5512)   build="T2Tv20"        ;;  # T2T v2.0 (CHM13 v1.1 w/ HG002 v2.7 Y) (chr name)
      84e78573982f3ea293bfeb54cd529309)   build="1K38p"         ;;  # Verily oddball GRCh38
   	  85c436650ffe85696c0fb51de4a3a74f)   build="THG1243v3"     ;;  # T2T HG01243 (aka PR1 Puerto Rican) (chr name)
   	  90814fe70fd8bbc59cacf2a3fd08e24c)   build="T2Tv09"        ;;  # T2T CHM13 v0.9
 	    a2fe6ab831d884104783f9be437ddbc0)   build="EBI38p"        ;;  # EBI GRCh38 Errored models
 	    a349c8b12bfffa22ea5326cfa878457f)   build="NCB37"         ;;  # Custom; identical to NCB37 but diff md5?
      a9634b94a29618dc3faf15a3060006ec)   build="HG18"          ;;  # hg 18 / NCBI 36
 	    b05113b52031beadfb6737bc1185960b)	  build="HG${mbuild}"   ;;  # HG 19 / NCBI 37
      b7884451f3069579e5f2e885582b9434)   build="1K38"          ;;  # 1K 38
 	    bbd2cf1448ccc0eaa2472408fa9d514a)   build="THGySeqp"      ;;  # ySeq HG38 w/ HG002 v2 Y
 	    bc811d53b8a6fc404d279ab951f2be4d)   build="HG17"          ;;  # hg 17 / NCBI 35
 	    bee8aebc6243ff5963c30abbd738d1f6)	  build="NCB38"         ;;  # NCBI 38 Genbank (orig top_level)
      c182b40ef3513ef9a1196881a4315392)   build="HPPv1"         ;;  # HPP CHM13 v1 with GRCh38 Y
 	    ca2e97bc5ecff43a27420eee237dbcc3)   build="EBI37p"        ;;  # EBI GRCh37 Errored models (extra Y)
      e9438f38ad1b9566c15c3c64a9419d9d)   build="T2Tv11"        ;;  # T2T CHM13 v1.1 (original)
      eec5eb2eeae44c48a31eb32647cd04f6)	  build="EBI38" 	      ;;  # EBI 38
      f7c76dbcf8cf8b41d2c1d05c1ed58a75)   build="NCB37"         ;;  # NCBI GRCh37 (RefSeq)
      d41d8cd98f00b204e9800998ecf8427e)   build="UNK"               # d41d8... is md5sum of empty file
        echo "$fileb: ***WARNING: Unknown Build (0 chromosomes)"
        errp+=" ***WARN: Uknown Build (0 chromosomes)***"                                          ;;
      *)                                  build="UNK"               # Anything else is unrecognized
        echo "$fileb: ***WARNING: Never encountered Build (${md5ff})"
        errp+=" ***WARN: Unrecognized Build Model***"                                          ;;
    esac
    build="${build}${msn}"			# Append sequence naming convention used

    # Create tab separated list of original primary SN names in fixed order (not upcased)
    # So cheat knowing only thing upcased is "chr" as all other primary name forms have upcase already
    chrSNt=$(cut -f1   "$filed.fptmp" | tr "\n" "\t" | tr "CHR" "chr" | sed 's/SN://g')

    # Save all the values as single line, tab separated text file; append to project file for directory
	  #  Note: chrSNt is already a multi-value, tab separated string that we simply use to create multiple columns
    printf $'\"%s\"\t%s\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\"%s\"\t%s\n' \
      "$filen" $build "$snct" "$md5bs" "$md5cs" "$md5fs" "$md5fp" "$md5sp" "$md5i" "$md5fc" \
      "$chrMSN" "$chrMLN" "$chrMM5" "$errp" "$chrSNt" > "$filed.wgse"

    # Final error check -- in case previous uncaught error left intermediate .wgse file or multi-line stats (mtDNA)
    wcnt=$(wc -l < "$filed.wgse")
    if (( wcnt != 1 )) ; then
      echo "$fileb: ***ERROR: Failed generating final WGSE stats file"
      $rmfx "$filed.wgse" || true
      exit_status+=1
    else
      cat "$filed.wgse" >> "${sumFPB}.csv"
    fi

    \rm "$filed.otmp" "$filed.stmp" "$filed.pstmp" "$filed.fptmp"  "$filed.mtmp" || true
  }
done

# Some basic stats for the Reference Genome study document when run on a directory with all known models (assume GNU utils)
if "$wholedir"; then
  $debugmode && echo "***DEBUG: Processing WGSE stats file for whole directory processing" && echo

  # To let us cd to dir then use local basename for wc; FPB always includes at least "./"; easier then sed on wc result)
  sumFP=$(dirname "$sumFPB")

  tail -n +2 "${sumFPB}.csv" > temp.csv     # Create Headless TSV file (strip first line of column labels)
  echo "Entries: " $(wc -l "temp.csv") | sed "s#temp#${sumFB}#"    # Use original file base name; not headless temp one
  cut -f2     "temp.csv" | sort -V | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_uniq_build.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_uniq_build.csv")
  cut -f8,13- "temp.csv" | sort --key=3,3 | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_uniq_ChrSN.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_uniq_ChrSN.csv")
  cut -f7     "temp.csv" | sort | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_uniq_ChrLNM5.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_uniq_ChrLNM5.csv")
  cut -f3     "temp.csv" | sort -n | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_uniq_SNcnt.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_uniq_SNcnt.csv")
  \rm temp.csv || true

  # For convenience, create a table of all Sequence entries in all FASTA's in the specified directory
  grep -v "^@HD" "$1"/*dict | cut -f2- | sort --key=2.4brn,2 -k3.4b,3 -k1.4b,1 > "${sumFPB}_dict.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_dict.csv")
  cut -f1-3 "${sumFPB}_dict.csv" | sort -V | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_dict_uniq_SNLNM5.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_dict_uniq_SNLNM5.csv")
  cut -f2-3 "${sumFPB}_dict.csv" | sort -rV | uniq -c | awk '{gsub(/^ +/,"") gsub(/ /,"\t")} {print $0}' > "${sumFPB}_dict_uniq_LNM5.csv"
  echo "Entries: " $(cd "${sumFP}" || true ; wc -l "${sumFB}_dict_uniq_LNM5.csv")

fi

(return 0 2>/dev/null) && return $exit_status || exit $exit_status
