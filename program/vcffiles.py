# coding: utf8
#
# VCF File Subsystem / module
#
# Part of the
# WGS Extract (https://wgse.bio/) system
#
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2023 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    VCFFile (class; vcffiles module) is the central handler for all things VCF.  Processing, determining stats,
    manipulating, etc.  The idea is a class instance for each VCF file specified.  Initially supporting multiple
    identified in the program (mostly as automatically found from the BAM specification).  But eventually
    raising VCFs to the same level and class as BAM (and eventually FASTQ also).  Once a VCF file(s) is selected
    in the interface or by the program automatically, stats on it will be collected and saved here for each file.

    VCF files are tougher in that there is even less information in the header about what they may contain. In lieu
    of scanning the whole file, we develop some quicker heuristics here as well (like done for BAM files).

    Although only created in Oct 2022, it contains replicated / modified code from the BAMfiles which has still some
    remnants from the original, single wgsextract.py file from v2.  So retained full copyright.
"""

import logging
import os       # for path, stat
from math import sqrt       # for process_bam_body
from vcf_parser import VCFParser

from utilities import  is_legal_path, nativeOS, universalOS, unquote, Error, Warning, wgse_message
from commandprocessor import run_bash_script, simple_command
from fastqfiles import determine_sequencer
import settings as wgse


######################################################################################################################
# BAM File Processing Error Classes
class VCFContentError(Error):
    def __init__(self, reason):
        self.reason = wgse.lang.i18n[reason]


class VCFContentErrorFile(Error):
    def __init__(self, reason, file):
        self.reason = wgse.lang.i18n[reason].replace('{{FILE}}', file)


class VCFContentWarning(Warning):
    def __init__(self, reason):
        self.reason = wgse.lang.i18n[reason]


######################################################################################################################
# VCF File Processing class

class VCFFile:
    """
        VCF File class object; including all stats and processing thereof
        In init, simply assert if error and let outer caller worry about deleting the object
        (i.e. we assume within invoked from within a try...except). Errors will be of class VCFContentError
        with an additional message parameter. The message parameter must be a LanguageStrings dictionary key.
    """

    def __init__(self, VCF_FN, nameContent):
        """
        True Class init that simply sets up the Class variables by initializing them (even if None).
        BAM file initialization is then called via _setup_BAM file at the end. Which can be called directly (internally)
        as well to setup a new BAM. But generally best to drop the Class instance and start over.

        nameContent is a subset of ({'SNP', 'InDel', 'CNV', 'SV'}, {'A', 'X', 'Y', 'M'})   # include "O"ther contigs?
        """
        self.file_oFN  = None   # Various renditions of the VCF file name, path, etc
        self.file_oFPB = None
        self.file_FN   = None
        self.file_qFN  = None
        self.file_FBS  = None
        self.file_FPB  = None
        self.file_FP   = None
        self.file_FB   = None
        self.file_FS   = None

        self.file_type = None   # for "VCF" or "BCF" (will we ever see BCF?)
        self.var_type  = set()  # set of {'SNP'. 'InDel', 'CNV', 'SV'} or less
        self.content = set()    # set of {'A', 'X', 'Y', 'M'} or less

        self.Header    = ""      # Will save complete VCF Header in text form (bcftools view -h)

        self.Refgenome = None    # {hg,GRCh}{19/37,38}, hs37d5, hs38DH (not all combinations valid)
        self.RefgenomeNew = None  # Following Reference Model study: hg, 1kg, ebi, ncbi
        self.Refgenome_qFN = None  # quoted File name of reference genome used for VCF

        self.RefMito   = None    # rCRS, Yoruba, RSRS
        self.SNTypeC   = None    # "Chr", "Num" or "Acc" (NCBI Accession #)
        self.SNTypeM   = None    # M or MT
        self.SNCount   = 0       # Count of '##CONTIG=<ID=' entries in Header
        self.Build     = None    # Major Build (19?, 36, 37, 38, 99 for T2T)

        self.Indexed   = False

        self._setup_VCF(VCF_FN, nameContent)

    def _setup_VCF(self, VCF_FN, nameContent):
        """ Check VCF File name / path is legal; fill in basic stats from VCF file if so. Assert on error. """
        VCF_oFN = nativeOS(VCF_FN)
        if not is_legal_path(VCF_oFN):
            raise VCFContentError('errVCFFileSpecialChars')
        if not os.path.isfile(VCF_oFN) or os.path.getsize(VCF_oFN) < 10000:    # file is bad
            raise VCFContentError('errVCFContent')

        # set_VCF_file_settings ... just doing file names. Most content set later while processing file.
        self.file_oFN  = VCF_oFN
        self.file_FN   = VCF_FN
        self.file_qFN  = f'"{self.file_FN}"'

        self.file_FP   = os.path.dirname(VCF_FN)
        self.file_FP  += '/' if self.file_FP[-1] != '/' else ''  # Assure trailing slash; universalOS so forward slash

        self.file_oFPB = os.path.splitext(VCF_oFN)[0]
        self.file_FPB  = universalOS(self.file_oFPB)

        self.file_FBS  = os.path.basename(VCF_FN)
        self.disp_FBS  = self.file_FBS if len(self.file_FBS) < 31 else f'{self.file_FBS[:12]} ... {self.file_FBS[-13:]}'

        self.file_FB, self.file_FS = self.file_FBS.split(os.extsep, 1)
        self.file_type = self.file_FS.upper.replace(".GZ", "")   # strip trailing .gz, if it exsits, and upcase
        if self.file_type not in ["VCF", "BCF"]:    # Internal call may pass an incorrect file
            raise VCFContentError('errVCFFileExtension')

        self.file_stats = os.stat(self.file_oFN)
        if self.file_stats.st_size < 10000:
            raise VCFContentError('errVCFFileEmpty')
        self.relfsize = max(0.01, self.file_stats.st_size / (45 * 10.0 ** 9))    # Relative size to 45GB (time adjust)

        if nameContent and len(nameContent) > 0 :
            self.var_type  = nameContent[0]  # set of {'SNP'. 'InDel', 'CNV', 'SV'} or less
            if len(nameContent) > 1:
                self.content = nameContent[1]

        self.process_vcf_header()       # Quick but may have error exceptions
        self.process_vcf_body()         # Little longer (10 seconds max) but key ones like read length, etc.

        wgse.VCFs += (self)             # May need within the below calls, so setup global setting now

    def process_vcf_header(self):
        """ Reads and stores VCF header. Sets SNCnt, and tries to determine reference model. """

        # Rewritten to process more of header; and process header in Python and not BASH AWK/GREP scripts
        bcftools = nativeOS(unquote(wgse.bcftoolsx_qFN))
        vcffile  = self.file_oFN

        self.Header = simple_command((bcftools, "view", "-H", "--no-PG", vcffile))

        if self.Header < 1000 or "##fileformat" not in self.Header:     # Command success check
            raise VCFContentError('errVCFHeader')

        for line in self.header:
            if line.startswith("##reference="):
                self.Reference = line

        # Determine SN names as given as opposed to relying on reference model to determine
        accession = ["##contig=<ID=CM0", "##contig=<ID=CP", "##contig=<ID=J0", "##contig=<ID=NC_"]
        self.SNTypeC = "Chr" if "##contig=<ID=chr" in self.Header else \
                       "Num" if "##contig=<ID=1\t" in self.Header else \
                       "Acc" if any(x in self.Header for x in accession) else \
                       "Unk"
        self.SNTypeM = "MT" if "##contig=<ID=MT" in self.Header or "##contig=<ID=chrMT" in self.Header else "M"
        self.SNCount = self.Header.count('##contig=<ID=')   # Almost sufficient for an analysis model determination
        logging.debug(f"VCF ID Type: {self.SNTypeC}, {self.SNTypeM}; VCF ID#:{self.SNCount}")
        
        self.indexed = self.check_for_vcf_index()

        self.determine_reference_genome()

    def process_vcf_body(self, reentrant=None):
        """
        Processes the BAM body.  Operate on first 100K sequence entries (2 million if Nanopore).
        Called internally from the "Show Stats" button as well.
        Will not have run idxstats yet if detected unsorted, unindexed or a CRAM.

        Used mainly to look at body of VCF to maybe determine content and type (to override initial setting
        determined from file name).
        """

    def determine_reference_genome(self):
        """
        Try to determine the reference genome from just the header information. VCF's have less than BAM's usually.

        The VCF Header is set from the caller creating the VCF.  It is generally carried through in full when subsets
        of the VCF are made (e.g. Y-only VCF).  Header "##contig=<ID=,length=" records should be based only on the
        reference model and not actual content of variants called. Dante includes the ##reference= entry but
        Nebula does not.

        Generic, major model naming comes from:
         hg if old style "chr22/chrM" sequence names, GRCh if numeric-only
         19 if old style hg with Yoruba, 37 if build 37 model lengths, 38 if Build 38 model lengths, ditto T2T 99
         Number of "contig" entries in the header is a good indicator of class and other model characteristics

        NOTE: Much of this code is duplicated in the BAMFiles class call of the same name.  Should likely pull out
        just the differences and make a single, generic call that adopts to the differences.  Maybe simply extract
        the SN/ID and LN/length fields into separate entry and process that?
        """
        from mainwindow import ask_refgenome

        self.Build, self.RefMito = \
            (99, "Custom") if any(x in self.Header for x in ["length=248387561", "length=248387328", "length=248387497",
                    "length=248415701", "length=62456832", "length=62460029", "length=62480187", "length=154343774",
                    "length=25843790",  "length=154259566", "length=154259625", "length=154269076", "length=154349815",
                    "length=154434329"]) else \
            (38, "rCRS")   if any(x in self.Header for x in ["length=248956422", "length=57227415", "length=156040895"]) else \
            (37, "rCRS")   if any(x in self.Header for x in ["length=249250621", "length=59373566", "length=155270560"]) else \
            (18, "Yoruba") if any(x in self.Header for x in ["length=247249719", "length=57772954", "length=154913754"]) else \
            (17, "Yoruba") if any(x in self.Header for x in ["length=245522847", "length=57701691", "length=154824264"]) else \
            (16, "Yoruba") if any(x in self.Header for x in ["length=246127941", "length=50286555", "length=153692391"]) else \
            (15, "Yoruba") if any(x in self.Header for x in ["length=245203898", "length=50961097", "length=152634166"]) else \
            ( 0, "Unknown")

        if self.Build == 37 and any(x in self.Header for x in ["M,length=16571", "MT,length=16571"]):
            self.Build, self.RefMito = (19, "Yoruba")
        # Todo handling RSRS model -- look for spacers?

        logging.debug(f"VCF Build: {self.Build:3d}")
        logging.debug(f"VCF Ref Genome Mito: {self.RefMito}")

        # New form is based on Major / Class mechanism in Reference Genome study: https://bit.ly/34CO0vj
        #  Used to rely on SNCount.  But oddball, unrecognized ref genomes in BAMs may have something close set by user.
        #  So, rely more on key items found in special ref genomes. Keep patching this.  Need to use the MD5sum method.
        if self.Build == 99:        # Used chr1, chrX or chrY lengths to determine Build99 earlier
            # Use X and Y length to determine T2T / HPP model
            if "length=62456832" in self.Header and "length=154343774" in self.Header:  # CHM13 v1.1 & HG002 v2 XY
                self.Refgenome = self.RefgenomeNew = "THGv20"
            elif "length=62460029" in self.Header and "length=154349815" in self.Header:  # CHM13 v1.1 & HG002 v2.7 XY
                self.Refgenome = self.RefgenomeNew = "THGv27"
            elif "length=62480187" in self.Header and "length=154434329" in self.Header:  # HG01243 v3 PR1 Puerto Rican
                self.Refgenome = self.RefgenomeNew = "THG1243v3"
            elif "length=57277415" in self.Header and "length=154259566" in self.Header:  # CHM13 v1.1 & GRCh38 Y
                self.Refgenome = self.RefgenomeNew = "HPPv11"
            elif "length=57277415" in self.Header and "length=154259625" in self.Header:  # CHM13 v1 X & GRCh38 Y
                self.Refgenome = self.RefgenomeNew = "HPPv1"
            elif "length=62456832" in self.Header and "length=156040895" in self.Header:  # GRCh38 w/ HG002 v2 Y  (ySeq)
                self.Refgenome = self.RefgenomeNew = "THGySeqp"
            elif "length=62460029" in self.Header and "length=154259566" in self.Header:  # CHM13 v1.1 & HG002 v2.7 Y
                self.Refgenome = self.RefgenomeNew = "T2Tv20"
            elif "length=154259566" in self.Header:  # Must be plain CHM13 v1.1 as special Y not found
                self.Refgenome = self.RefgenomeNew = "T2Tv11"
            elif "length=154259625" in self.Header:  # Must be plain CHM13 v1 as special Y not found
                self.Refgenome = self.RefgenomeNew = "T2Tv10"
            elif "length=154259664" in self.Header:  # Must be plain CHM13 v0.9
                self.Refgenome = self.RefgenomeNew = "T2Tv09"
            else:  # Think it is T2T / HPP build but cannot recognize Y or X length
                logging.debug(f'Unrecognized T2T model; no matching X or Y length in BAM to known models.')
        elif "##contig=<ID=NC_007605" in self.Header and self.SNCount in [85, 86]:
            # hs37 class have SN:NC007605 (EBV in Numeric naming) sans human_g1k / GRCh37.primary_assembly models
            # human_1kg model is 84 and one less SN than hs37.fa.gz; hs37d5 is 86 and one more than hs37 (hs37d5 decoy)
            self.Refgenome = "hs37d5" if "##contig=<ID=hs37d5" in self.Header else "hs37"
            self.RefgenomeNew = "1k37g"
        elif "##contig=<ID=chrEBV" in self.Header and self.SNCount in [85, 298]:
            self.Refgenome = None     # 1KGenome analysis models in NCBI Genbank Archive or UCSC; but not handled yet
            self.RefgenomeNew = "1k37h"
        elif self.SNCount == 84:      # human_g1k (if Num), GRCh37.primary_assembly.genome (if Chr)
            self.Refgenome = "hs37-" if self.SNTypeC == "Num" else \
                             "GRCh37-" if self.SNTypeC == "Chr" else None
            self.RefgenomeNew = "1K37g" if self.SNTypeC == "Num" else \
                                "EBI37h" if self.SNTypeC == "Chr" else None
        elif ("##contig=<ID=chrEBV" in self.Header or "##contig=<ID=EBV" in self.Header) and \
              self.SNCount in [195, 456, 2581, 3366]:
            # hs38DH is 3366 SN count and has SN:HLA- unique; hs38 is 195 and hs38a is 456; all uniquely have chrEBV
            # hs38s is sequencing.com model of GCA_000001405.15_GRCh38_no_alt_plus_hs38d1.fna.gz with 22_KI270879v1_alt
            self.Refgenome = "hs38DH"  if "##contig=<ID=chr22_KI270879v1_alt" in self.Header and self.SNCount == 3366 else \
                             "hs38s"   if "##contig=<ID=22_KI270879v1_alt" in self.Header and self.SNCount == 2581 else \
                             "hs38a"   if self.SNCount == 456 else \
                             "hs38"  # if self.SNCount == 195
            self.RefgenomeNew = "1k38" + ("pg" if self.Refgenome == "hs38s" else "h")
        elif self.SNCount in [25, 93, 297, 455, 639]:    # SNCounts of 6 hg/ebi models; 2 duplicated
            # Handle the hg (UCSC) and GRCh (EBI) models here
            self.Refgenome  = "hg" if self.SNTypeC == "Chr" else "GRCh"
            self.Refgenome += str(self.Build)   # Build already takes into account mito for 19/37 differentiation
            self.RefgenomeNew = self.Refgenome.replace("GRCh", "EBI") + ("g" if self.SNTypeC == "Chr" else "h")
        logging.debug(f'Ref Genome: {self.Refgenome}, by New nomenclature: {self.RefgenomeNew}')

        # We give up if still not set; simply ask the user
        if not self.Refgenome:  # not set yet: and (self.chrom_types["A"] > 1 or self.chrom_types["Y"] > 1):
            ask_refgenome(self, type="VCF", inBAM=True)  # Return value in self class pointer
            logging.debug(f'Ref Genome (User): {self.Refgenome}')
            # todo ask user to send VCF header so we can get its signature for future runs

        self.Refgenome_qFN = wgse.reflib.get_refgenome_qFN(self.Refgenome)
        tempFN = unquote(self.Refgenome_qFN)
        logging.debug(f'Ref Genome File: {tempFN}')

    def chrom_types_str(self):
        """
        Create a string of comma-separated BAM SN components: Auto(somes), X, Y, Mito, etc of BAM / CRAM for stats
        Typical 30x WGS male BAM file would return: "Auto, X, Y, Mito, Unmap, Other" (in English lang setting)
        """
        count = 0
        chrom_types_str = ""
        for key, value in self.chrom_types.items():
            if value > 1 and not (key == 'Y' and self.gender == 'Female'):  # Do not show Y on Female
                count += 1
                chrom_types_str += (", " if chrom_types_str else "") + wgse.lang.i18n[key]
                ''' Key values are one of: ("A", "X", "Y", "M", "O"). Single character lookup in language file. '''
        if count == 1:
            chrom_types_str += f' {wgse.lang.i18n["Only"]}'
        return chrom_types_str

    def refgenome_str(self):
        """ Create string of ref genome, mito genome, and SN Count of VCF for stats """

        if self.Refgenome:
            result = self.Refgenome
        else:
            bldstr = f' {str(self.Build)}' if self.Build else ""
            result = wgse.lang.i18n["Unknown"] + bldstr
        if self.SNTypeC:
            result += f' ({self.SNTypeC})'
        if self.RefMito:
            result += ", " if result else ""
            result += self.RefMito
        if self.SNCount:
            result += ", " if result else ""
            result += f'{self.SNCount} {wgse.lang.i18n["SNs"]}'
        return result

    def filestatus_str(self):
        """ Create string of sorted, indexed and file size status of BAM / CRAM for stats """
        result = ""
        #if self.Indexed:
        #  result += ", " if result else ""
        result += wgse.lang.i18n["Indexed"] if self.Indexed else wgse.lang.i18n["Unindexed"]
        # File size
        result += ", " if result else ""
        result += f'{round(self.file_stats.st_size / (10.0 ** 9), 1)} {wgse.lang.i18n["GBs"]}' if self.file_stats.st_size >= 10 ** 9 else \
                  f'{round(self.file_stats.st_size / (10.0 ** 6), 1)} {wgse.lang.i18n["MBs"]}' if self.file_stats.st_size >= 10 ** 6 else \
                  f'{round(self.file_stats.st_size / (10.0 ** 3), 1)} {wgse.lang.i18n["KBs"]}' if self.file_stats.st_size >= 10 ** 3 else \
                  f'{round(self.file_stats.st_size, 1)} {wgse.lang.i18n["Bs"]}'
        return result

    def get_chr_name(self,chrtype):
        """ Return chromosome name of "type" (MT or Y) specific to the BAM reference model. """
        if self.SNTypeC == "Chr":
            chrM = "chrM" if self.SNTypeM == "M" else "chrMT"
            chrY = "chrY"
        elif self.SNTypeC == "Num":
            chrM = "MT" if self.SNTypeM == "MT" else "M"
            chrY = "Y"
        elif self.SNTypeC == "Acc":
            chrM = "M"
            chrY = "Y"
        else:
            chrM = chrY = ""
        return chrM if chrtype == "M" else chrY

    def check_for_vcf_index(self):
        """ VCF Index file exists check.  """
        # todo need to add check in OUT directory, if set. And add option to use in samtools calls if located there
        self.Indexed = os.path.isfile(self.file_oFN + ".tbi") or os.path.isfile(self.file_oFN + ".TBI")
        # todo case sensitive check needed? how does bcftools handle?
        return self.Indexed

'''  Part of BAMfiles (copied here when file duplicated)
    def find_VCFs(self):
        """
        A special to search for commonly-named VCFs from a vendor for the loaded BAM file.  As a precursor to creating
        them from the BAM.  In here as if they are found, will set class variable so they can be used and associated.
        Eventually, add "get_" routine like for reference_genome to ask the user if not found. Saves time and disk
        space to reuse instead of recreating. Helpful especially when we start saving and restore BAM/FASTQ/VCF
        library stats for a complete test kit.
        """
        # todo checking dates (compared to BAM) and sizes (minimum)

        # First check if names have already been set in this BAM class
        if not (self.VCFs and len(self.VCFs) > 0):
            self.VCFs = VCFs = {}

            # First look in output directory for VCFs we create (our naming convention)
            if wgse.outdir and wgse.outdir.FPB:
                VCFs += {
                    f'{wgse.outdir.FPB}_SNP.vcf.gz'  : ({'SNP'}, {'A', 'X', 'Y', 'M'}),
                    f'{wgse.outdir.FPB}_InDel.vcf.gz': ({'InDel'}, {'A', 'X', 'Y', 'M'}),
                    f'{wgse.outdir.FPB}_CNV.vcf.gz'  : ({'CNV'}, {'A', 'X', 'Y', 'M'}),
                    f'{wgse.outdir.FPB}_SV.vcf.gz'   : ({'SV'}, {'A', 'X', 'Y', 'M'}),
                    f'{wgse.outdir.FPB}.vcf.gz'      : ({'unk'}, {'A', 'X', 'Y', 'M'})
                }

            # Now look in BAM directory for common test company names that come with their BAMs
            VCFs += {
                f'{self.file_FPB}.filtered.snp.vcf.gz'    : ({'SNP'}, {'A', 'X', 'Y', 'M'}),    # Dante Labs
                f'{self.file_FPB}.filtered.indel.vcf.gz'  : ({'InDel'}, {'A', 'X', 'Y', 'M'}),  # Dante Labs
                f'{self.file_FPB}.cnv.vcf.gz'             : ({'CNV'}, {'A', 'X', 'Y', 'M'}),  # Dante Labs & Sequencing
                f'{self.file_FPB}.sv.vcf.gz'              : ({'SV'}, {'A', 'X', 'Y', 'M'}),   # Dante Labs & Sequencing
                f'{self.file_FPB}.vcf.gz'                 : ({'SNP', 'InDel'}, {'A', 'X', 'Y', 'M'}),  # Nebula & ySeq
                f'chrY_called_{self.file_FPB}.vcf.gz'     : ({'SNP'}, {'Y',}),                    # ySeq
                f'chrY_cleaned_{self.file_FPB}.vcf.gz'    : ({'SNP'}, {'Y',}),                    # ySeq
                f'chrY_derived_{self.file_FPB}.vcf.gz'    : ({'SNP'}, {'Y',}),                    # ySeq
                f'chrY_INDELS_{self.file_FPB}.vcf.gz'     : ({'InDel'}, {'Y',}),                  # ySeq
                f'chrY_novel_SNPs_{self.file_FPB}.vcf.gz' : ({'SNP'}, {'Y',}),                    # ySeq
                f'{self.file_FPB}.snp-indel.genome.vcf.gz': ({'SNP', 'InDel'}, {'A', 'X', 'Y', 'M'}),  # Sequencing.com
                f'{self.file_FPB}.mito.vcf.gz'            : ({'SNP', 'InDel'}, {'A', 'X', 'Y', 'M'})   # Sequencing.com
                # Full Genomes Corp?
            }
            for key, val in VCFs.items():
                if os.path.exists(nativeOS(key)):
                    self.VCFs += { key, val }

        if not (self.VCFs and len(self.VCFs) > 0):
            # todo ask the user if they know where they are and let them set.  Implement once VCF stats save happens.
            pass

        process_vcf_headers(self)

        return self.VCFs
'''
