#!/bin/bash
# run-Gen-Ref script originally from Heng Li's BWA-KIT release to make the 5 1K Genome project reference genomes ("hs")
# See https://github.com/lh3/bwa/tree/master/bwakit
# Updated for current paths and to optimize downloads; added all option
# Updated to utilize files of the bwa-kit/resource-GRCh38 folder from $reflibdir.
# Note: does not index and uses gzip instead of bgzip. But process_refgenomes.sh fixes all that.

# todo add standard WGSE zcommon.sh source directive here
root=${reflibdir}

url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
url37d5="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

if [ $# -eq 0 ]; then
	echo "Usage: $0 <hs38|hs38a|hs38DH|hs37|hs37d5|all>"
	echo "Analysis sets:"
	echo "  hs38     primary assembly of GRCh38 (incl. chromosomes, unplaced and unlocalized contigs) and EBV"
	echo "  hs38a    hs38 plus ALT contigs"
	echo "  hs38DH   hs38a plus decoy contigs and HLA genes (recommended for GRCh38 mapping)"
	echo "  hs37     primary assembly of GRCh37 (used by 1000g phase 1) plus the EBV genome"
	echo "  hs37d5   hs37 plus decoy contigs (used by 1000g phase 3)"
	echo ""
	echo "Note: This script downloads human reference genomes. For hs38a and hs38DH, it needs additional"
	echo "      sequences and ALT-to-REF mapping included in the bwa.kit package."
	exit 1;
fi

if [ "$1" == "hs38DH" ]; then
	(wget -O- $url38 | gzip -dc; cat "$root"/hs38DH-extra.fa) > "$1".fa
	[ ! -f "$1".fa.alt ] && cp "$root"/hs38DH.fa.alt "$1".fa.alt
elif [ "$1" == "hs38a" ]; then
	wget -O- $url38 | gzip -dc > "$1".fa
	[ ! -f "$1".fa.alt ] && grep _alt "$root"/hs38DH.fa.alt > "$1".fa.alt
elif [ "$1" == "hs38" ]; then
	wget -O- $url38 | gzip -dc | awk '/^>/{f=/_alt/?0:1}f' > "$1".fa
elif [ "$1" == "hs37d5" ]; then
	wget -O- $url37d5 | gzip -dc > "$1".fa 2>/dev/null
elif [ "$1" == "hs37" ]; then
	wget -O- $url37d5 | gzip -dc 2>/dev/null | awk '/^>/{f=/>hs37d5/?0:1}f' > "$1".fa
elif [ "$1" == "all" ]; then
	(wget -O- $url38 | gzip -dc; cat "$root"/hs38DH-extra.fa) > hs38DH.fa
	[ ! -f hs38DH.fa.alt ] && cp "$root"/hs38DH.fa.alt hs38DH.fa.alt
	wget -O- $url38 | gzip -dc > hs38a.fa
	[ ! -f hs38a.fa.alt ] && grep _alt "$root"/hs38DH.fa.alt > hs38a.fa.alt
	wget -O- $url38 | gzip -dc | awk '/^>/{f=/_alt/?0:1}f' > hs38.fa
	wget -O- $url37d5 | gzip -dc > hs37d5.fa 2>/dev/null
	wget -O- $url37d5 | gzip -dc 2>/dev/null | awk '/^>/{f=/>hs37d5/?0:1}f' > hs37.fa
else
	echo "ERROR: unknown genome build"
fi

[ ! -f "$1".fa.bwt ] && echo -e "\nPlease run 'bwa index $1.fa'...\n"

