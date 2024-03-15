# coding: utf8
# Copyright (C) 2020-2022 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    ReferenceLibrary module in progress.  Eventually handling the dynamic download of any needed reference library
    element (genome reference model, GFF/GTF variant models, liftover file, etc). Will always require a seed of the
    needed files (where to find them, characteristics to select them, etc)  Originally setup as a standalone program.

    Started as the accurate identification of the  reference genome used to create the BAM.  Betav1/2 of the tool
    used a simple determination of the main, likely reference genome that was actually inaccurate for the models
    included with it. v3 now expands to more exactly identify the model used by using a hash signature created from
    the model FA files and matched to a similar signature derived from the BAM file header. This is not the same as
    the md5 signature used by SAMTools for each sequence name.  Instead, it is am md5 signature of the BAM header
    SQ entries (SN's and lengths found in a BAM header and created from the FA's .dict).  See https://bit.ly/34CO0vj
    for the companion document developed with the work here.
"""
import os.path
# from tkinter.ttk import Button, Label
from tkinter import Toplevel, Radiobutton, Label, StringVar
try:
    from tkmacosx import Button
except ImportError:
    from tkinter import Button

from utilities import DEBUG, nativeOS, is_legal_path, unquote, wgse_message
from commandprocessor import run_bash_script
import settings as wgse


class ReferenceLibrary:
    """
    In preparation for a full function Reference Library capability (library maintenance of multiple types of objects;
    manipulation, automatic keep-up-to-date, let users add entries, etc), creating this class now to start.
    For now, mostly to capture and handle the Reference Library file-name / access variables that were in settings.
    Too much functionality was buried in mainwindow that needed to be pulled out to handle the settings storage and
    restore.  So better to put it into a class here than mainwindow which should be UI focused.
    """
    
    # Todo Liftover files: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ and
    #  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ should be added to reference library
    #  as should VCF annotations, BED files, etc

    # Todo read in genomes.csv to create this table and more
    _RefGenome_FBS = {
        'hg19': "hg19.fa.gz",
        'hg37': "hg19_yseq.fa.gz",
        'hg37w': "hg19_WGSE.fa.gz",
        'hs37-': "human_g1k_v37.fasta.gz",
        'hs37': "hs37.fa.gz",
        'hs37d5': "hs37d5.fa.gz",
        'GRCh37-': "GRCh37.primary_assembly.genome.fa.gz",
        'GRCh37': "Homo_sapiens.GRCh37.dna.toplevel.fa.gz",
        "hg38": "hg38.fa.gz",
        "hs38": "hs38.fa.gz",
        "hs38a": "hs38a.fna.gz",
        "hs38d1": "hs38d1.fna.gz",
        "hs38d1s": "hs38d1s.fa.gz",
        "hs38d1v": "GRCh38_Verily_v1.genome.fa.gz",
        "hs38d1a": "hs38d1a.fna.gz",
        "hs38DH": "hs38DH.fa.gz",
        "GRCh38": "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
        "GRCh38-": "GRCh38.primary_assembly.genome.fa.gz",
        "THGv27": "hg002xy_v2.7.fasta.gz",
        "THGv20": "hg002xy_v2.fasta.gz",
        "HPPv11": "CHM13v11Y.fa.gz",
        "HPPv1": "CHM13v1Y.fa.gz",
        "THG1243v3": "hg01243_v3.fna.gz",
        "THGySeqp": "hg38_CP086569.fasta.gz",
        "T2Tv20": "chm13v2.0.fa.gz",
        "T2Tv20a": "GCA_009914755.4.fna.gz",
        "T2Tv11": "chm13.draft_v1.1.fasta.gz",
        "T2Tv10": "chm13.draft_v1.0.fasta.gz",
        "T2Tv09": "chm13.draft_v0.9.fasta.gz",
        "unknown": "unknown"
    }

    def __init__(self, default_FP):
        self.valid = False      # If current reference library setting is valid (False if default and does not exist)
        self.set = None         # If set by settings file or user; default has been overriden
        self.default_FP = None  # Default Reference Library location FP
        self.FP  = None         # Reference Library File (Directory) Pointer
        self.oFP = None
        self.gen_FP  = None     # (Reference) Genome File (Directory) Pointer
        self.gen_oFP = None
        self.FB  = None         # Final (leaf) directory name of Reference Library location (e.g. "reference")
        self.cma_FP = None      # Chip MicroArray Templates and support File (Directory) Pointer
        self.cma_oFP = None
        self.liftover_chain_oFN = None      # Liftover Chain pointer (currently Build38 -> Build37 only)
        self.askRefgenWindow = None

        self.change(default_FP)     # Call Reference Library location Change with Default location (for initial setup)
        # if not self.FP:           # Default will not exist if installation made with reflib already redirected
        #     Setting Reference Library default failed; really fatal and should probably exit (raise exception)
        #     DEBUG(f'***ERROR: Cannot set default Reference Library directory: {default_FP}')
        #     del self

    def change(self, dir_FP):

        if not dir_FP or dir_FP == "/" or dir_FP == "C:/" or dir_FP == "/cygdrive/c/":
            # Directory must not be the root file system (also check for null value early on)
            if wgse.lang.i18n:         # Only report full issue if wgse.lang setup
                wgse_message("error", 'InvalidRefLibTitle', True,
                             wgse.lang.i18n['errRefLibPath'].replace("{{DIR}}", dir_FP))
            else:
                DEBUG(f'***WARNING: Empty path or Root directory Reference Library not allowed: "{dir_FP}"')
            return
        elif dir_FP == self.FP:  # No change. Silently return. If initilization, self.FP not set yet
            return

        dir_FP += '/' if dir_FP[-1] != '/' else ''  # Assure trailing slash; universalOS so always forward slash
        dir_oFP = nativeOS(dir_FP)

        gen_FP  = f'{dir_FP}genomes/'
        gen_oFP = nativeOS(gen_FP)
        cma_FP = f'{dir_FP}microarray/'
        cma_oFP = nativeOS(cma_FP)

        # Error check before setting class instance values (language not setup if default reference library)
        self.valid = True
        if not(is_legal_path(dir_oFP) and os.path.isdir(dir_oFP)):
            # Simply a bad path (not a directory, does not exist, etc)
            self.valid = False
            if wgse.lang.i18n and self.default_FP:
                wgse_message("error", 'InvalidRefLibTitle', True,
                             wgse.lang.i18n['errRefLibPath'].replace("{{DIR}}", dir_oFP))
                return
            else:
                # If WGS Extract installed with a .wgsextract reflib setting already redirecting the reflib, then
                #   the default installation reference/ directory will not exist.  So just issue a warning but
                #   still proceed as if valid so default is set and GUI operates accordingly.
                DEBUG(f'***WARNING: Reference Library (default) directory does not exist: "{dir_FP}"')
        elif not(os.path.isdir(gen_oFP) and os.path.isdir(cma_oFP)):
            # We hope they have setup a library with content; at minimum, we need our special directories there
            self.valid = False
            if wgse.lang.i18n and self.default_FP:
                wgse_message("warning", 'InvalidRefLibTitle', False, 'errRefLibEmpty')
            else:
                DEBUG(f'***WARNING: Reference Library (default) directory exists but is missing content: "{dir_FP}"')
            # return   # No return; just warning and continue on

        # Now that we have error checked (confirmed) values, go ahead and store in the Class instance
        if not self.default_FP:                     # If first time initializsing, set the default
            self.default_FP = dir_FP
        self.set = dir_FP != self.default_FP        # Meaning set by user or settings; not default
        self.FP  = dir_FP       ; self.oFP = dir_oFP
        self.gen_FP  = gen_FP   ; self.gen_oFP = gen_oFP
        self.cma_FP = cma_FP    ; self.cma_oFP = cma_oFP
        self.FB = os.path.basename(dir_FP[:-1])  # Trick; remove trailing slash so basename returns directory name

        self.liftover_chain_oFN = f'{self.oFP}hg38ToHg19.over.chain.gz'
        # Todo have chain files for T2T; can we make microarray work by creating T2T chain files then liftover result?

        # self.snlookup =

        DEBUG(f'New Reflib Directory: {self.FP}')

    def get_refgenome_qFN(self, ref_genome):
        """
        This simplistic extension of the original does not recognize NCBI models nor minor patch variants
        It does continue to include the original models (even if in error and should not be used) to aid
        in conversion back to something useful.
        """

        RefGenomeFBS = self._RefGenome_FBS.get(ref_genome, "unknown")
        if ref_genome and "unknown" not in [ref_genome, RefGenomeFBS]:
            return f'"{self.gen_FP}{RefGenomeFBS}"'
        # else return null

    def missing_refgenome(self, refgenome_qFN, parent=wgse.window):
        """ Check stored reference genome filename for existing and large enough; report error and give chance to
        correct if not; return final status. Calling this routine assumes the ref genome is required.  """
        missing = True      # Default return value

        if not refgenome_qFN:   # Nothing specified; maybe user cancelled specifying; simply return
            return missing

        refgenome_FN  = unquote(refgenome_qFN)
        refgenome_oFN = nativeOS(refgenome_FN)
        refgenome_FBS = os.path.basename(refgenome_FN)

        refgenome = [key for key in self._RefGenome_FBS if self._RefGenome_FBS[key] == refgenome_FBS][0]

        # Also check for the post-process file from our application
        special = refgenome_oFN.replace("fa.gz", "wgse").replace("fna.gz", "wgse").replace("fasta.gz", "wgse")

        # Standard case -- refgenome file found and correct size
        if os.path.isfile(refgenome_oFN) and os.path.getsize(refgenome_oFN) > 500000000 and \
           os.path.isfile(special):
            return not missing

        # Give user a chance to resolve issue by fetching the missing file; check again after
        message = wgse.lang.i18n["errRefGenFileMissing"]
        message = message.replace("{{FILE}}", refgenome_FBS).replace("{{REFGEN}}", refgenome)
        result = wgse_message("yesnocancel", "errRefGenFileMissingTitle", True, message)

        if result is None:  # User hit Cancel to decline to fix the missing refgenome. Simply return with missing.
            return missing

        if result:          # User hit Yes to have us to load the missing reference genome for them
            get_and_process = nativeOS(f'{wgse.install_FP}scripts/get_and_process_refgenome.sh')
            command = [wgse.bashx_oFN, get_and_process, refgenome_FBS, wgse.prefserver]
            run_bash_script("GenLoadGenome", command, parent=parent, direct=True)
            # The above generates a one line command to execute directly (no temp folder bash script)
        # else not result: user hit No after running Library.* themself; so hopefully should be available now ...

        # We need to check again in case the refgenome now exists (either we loaded or user did)
        if os.path.isfile(refgenome_oFN) and os.path.getsize(refgenome_oFN) > 500000000 and \
           os.path.isfile(special):
            return not missing

        # Simply do an error pop-up to let them know the current command will cancel (as the file is still missing)
        message = wgse.lang.i18n["errRefGenFileStillMissing"]
        message = message.replace("{{FILE}}", refgenome_FBS).replace("{{REFGEN}}", refgenome)
        wgse_message("error", "errRefGenFileMissingTitle", True, message)

        return missing

    @staticmethod
    def refgen_liftover(ref_genome):
        """
        Return companion reference genome name and SN Count from the reference genome name (only) passed in.
        """
        # We do not yet support NCBI models directly (the 4th class); we realign draft, early T2T style models to its final
        return  ("hg37",    93) if ref_genome == "hg38" else \
                ("hs37d5",  86) if ref_genome in ["hs38", "hs38a", "hs38d1", "hs38d1a", "hs38d1s","hs38DH", "hs38d1v"] else \
                ("GRCh37", 297) if ref_genome in ["GRCh38", "GRCh38-"] else \
                ("hg38",   455) if ref_genome in ["hg37", "hg19"] else \
                ("hs38",   195) if ref_genome in ["hs37d5", "hs37", "hs37-", "T2Tv20"] else \
                ("GRCh38", 639) if ref_genome in ["GRCh37", "GRCh37="] else \
                ("T2Tv20",  25) if ref_genome in ["T2Tv11", "T2Tv10", "T2Tv09", "THGySeqp", "THGv27", "THGv20",
                                                  "THG1243v3", "HPPv11", "HPPv1"] else \
                ("error", 0)

    def get_reference_vcf_qFN(self, build, sqname, type="FullGenome"):
        if type == "Yonly":  # Y-only Variants from yBrowse DB; updated daily.  Thomas Krahn suggest liftover to Build38
            FBS = "snps_hg38.vcf.gz"   if  build == 38 and sqname == "Chr" else \
                  "snps_hg19.vcf.gz"   if (build == 37 or build == 19) and sqname == "Chr" else \
                  "snps_grch38.vcf.gz" if  build == 38 and sqname == "Num" else \
                  "snps_grch37.vcf.gz" if (build == 37 or build == 19) and sqname == "Num" else \
                  "error"  # Unknown
            return f'"{self.FP}{FBS}"'  # In reference_library; from ybrowse.org
        elif type == "Microarray":  # Special list of CombinedKit Variants (union of all generated kits)
            FBS = "All_SNPs_hg38_ref.tab.gz"   if  build == 38 and sqname == "Chr" else \
                  "All_SNPs_hg19_ref.tab.gz"   if (build == 37 or  build == 19) and sqname == "Chr" else \
                  "All_SNPs_GRCh38_ref.tab.gz" if  build == 38 and sqname == "Num" else \
                  "All_SNPs_GRCh37_ref.tab.gz" if (build == 37 or  build == 19) and sqname == "Num" else \
                  "error"  # Unknown
            return f'"{self.cma_FP}{FBS}"'  # Found in the Microarray folder
        else:  # type == "FullGenome" (use dbSNP)      # Found in reference_library; from https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
            pass  # not yet implemented nor called; these files are very large. 15GB each or larger. should only be used for annotation

    @staticmethod
    def get_wes_bed_file_qFN(build, chrtype, ymt_only=False):
        """
        These WES BED Files have been created / modified from the two originals.
        (1) The numeric naming versions have been created.
        (2) The Mitochondria entries have been moved to the end to match chromosome order in the 1KGenome ref models
        (3) The xgen Build 38 files have had the mtDNA added (took same regions from TruSeq as all rCRS)
        (4) The bed files have been sorted to match Dante/Nebula sort order (1-22,X,Y,M) (sort -V with MT move)
        These two WES BED Files were chosen as they are based on the 1DT / TruSeq library prep and flow supported
        by Illumina and used by Dante.  The xgen file being the Build 38 version of the TruSeq Build 37 one.  Oddly,
        1/4 of the Y chromosome appears as zero coverage and so throws off stats.  Is it an N region covered by the
        BED file?  Analysis models have it nulled out for BWA alignment.
        (5) Technically, hg19 is Yoruba.  But we have kept the same mtDNA regions.

        The original BED files are found at:
        https://support.illumina.com/downloads/truseq-exome-product-files.html
        https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=3803
        The Y-only BED file started from an hg38 spreadsheet of regions from David Vance.  Obtained at
        https://docs.google.com/file/d/1_YwP_GRvSjFIDw5B5mhi0FeRlc2z2nLu/edit?usp=docslist_api&filetype=msexcel
        """
        wes_bed_files = {
            "19Chr": "TruSeq_Exome_TargetedRegions_v1.2.bed",   # Todo need to create a unique hg19 with Yoruba BED
            "37Chr": "TruSeq_Exome_TargetedRegions_v1.2.bed",
            "37Num": "TruSeq_Exome_TargetedRegions_v1.2num.bed",
            "38Chr": "xgen_plus_spikein.GRCh38.bed",
            "38Num": "xgen_plus_spikein.GRCh38num.bed"
        }
        ymt_only_bed_files = {
            "19Chr": "CombBED_McDonald_Poznik_Merged_hg37.bed", # Todo need to create a unique hg19 with Yoruba BED
            "37Chr": "CombBED_McDonald_Poznik_Merged_hg37.bed",
            "37Num": "CombBED_McDonald_Poznik_Merged_hg37num.bed",
            "38Chr": "CombBED_McDonald_Poznik_Merged_hg38.bed",
            "38Num": "CombBED_McDonald_Poznik_Merged_hg38num.bed"
        }
        """ymt_only_bed_files = {
            "19Chr": "BigY3_hg37.bed", # Todo need to create a unique hg19 with Yoruba BED
            "37Chr": "BigY3_hg37.bed",
            "37Num": "BigY3_hg37num.bed",
            "38Chr": "BigY3_hg38.bed",
            "38Num": "BigY3_hg38num.bed"
        }"""
        if ymt_only:
            return f'"{wgse.reflib.FP}{ymt_only_bed_files.get(str(build) + chrtype, "unknown")}"'
        else:
            return f'"{wgse.reflib.FP}{wes_bed_files.get(str(build) + chrtype,"unknown")}"'

    def ask_reference_genome(self, inBAM=True):
        """
        If we cannot figure the reference genome from BAM, ask the user.  Gives them option to say "unknown".
        Also used when need to ask for a reference genome during general alignment or cannot find a match for realign.
        Returns a typle of Refgenome code / string and the SN count for that reference
        """
        # Todo rework to use genomes.csv file to create radio button table or drop down; disabling ones with SN count
        #  smaller in ref than BAM here; wrong SN naming convention, etc. Biggest issue is patched models. Only
        #  latest GRCh models from Ensembl here include patch releases of base models.

        font = wgse.fonts.table if wgse and wgse.fonts else {'14': ("Times New Roman", 14)}

        # mainWindow.withdraw()
        rbRefgenome = StringVar(wgse.window)

        # Technically, may not be mainwindow that we are forked from; could be anywhere. But have no current window.
        askRefgenWindow = Toplevel(wgse.window)
        askRefgenWindow.transient()
        askRefgenWindow.protocol("WM_DELETE_WINDOW", self.ask_reference_genome_done)
        askRefgenWindow.title(wgse.lang.i18n['SelectReferenceGenome'])
        # askRefgenWindow.geometry("")
        askRefgenWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
        askRefgenWindow.columnconfigure(0, weight=1)
        askRefgenWindow.rowconfigure(0, weight=1)

        textToDisplay = ""
        if inBAM:  # Failed trying to determine Reference Genome of BAM; so ask user
            textToDisplay = wgse.lang.i18n['CouldntDetermineReferenceGenome'].replace("{{RefGen}}", wgse.BAM.refgenome_str())
            textToDisplay += wgse.lang.i18n['CautionSelectTheCorrectRefGenome'] + '\n'

        textToDisplay += wgse.lang.i18n['PleaseSelectReferenceGenome']
        Label(askRefgenWindow, text=textToDisplay, font=font['14'],
              anchor="w", justify="left").grid(column=0, row=0, columnspan=2, padx=5, pady=3)
        rowcnt = 1

        # Setup series of radio buttons for available reference genomes; disable buttons that clearly cannot
        #  be (incorrect Build or Chr Name type)
        rbhs37d5 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS37D5'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("hs37d5", 86))
        rbhs37d5.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhs37m = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHumanG1k'], font=font['14'],
                              indicatoron=False, variable=rbRefgenome, value=("hs37-", 84))
        rbhs37m.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbhs38 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38'], font=font['14'],
                             indicatoron=False, variable=rbRefgenome, value=("hs38", 195))
        rbhs38.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhs38dh = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38DH'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("hs38DH", 3366))
        rbhs38dh.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbhg37 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG37'], font=font['14'],
                             indicatoron=False, variable=rbRefgenome, value=("hg37", 93))
        rbhg37.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhg19 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG19'], font=font['14'],
                             indicatoron=False, variable=rbRefgenome, value=("hg19", 93))
        rbhg19.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbhg38 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG38'], font=font['14'],
                             indicatoron=False, variable=rbRefgenome, value=("hg38", 455))
        rbhg38.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhg37wg = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG37wgse'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("hg37wg", 25))
        rbhg37wg.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbgrch37 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsGRCh37'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("GRCh37", 297))
        rbgrch37.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbgrch38 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsGRCh38'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("GRCh38", 639))
        rbgrch38.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbt2tv20 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsT2Tv20'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("T2Tv20", 25))
        rbt2tv20.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbthg1243v3 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsTHG1243v3'], font=font['14'],
                                  indicatoron=False, variable=rbRefgenome, value=("THG1243v3", 89))
        rbthg1243v3.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbthgv27 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsTHGv27'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("THGv27", 25))
        rbthgv27.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbthgv20 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsTHGv20'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("THGv20", 25))
        rbthgv20.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbhppv11 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHPPv11'], font=font['14'],
                               indicatoron=False, variable=rbRefgenome, value=("HPPv11", 24))
        rbhppv11.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhppv1 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHPPv1'], font=font['14'],
                              indicatoron=False, variable=rbRefgenome, value=("HPPv1", 24))
        rbhppv1.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        rbhs38d1 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38d1'], font=font['14'],
                            indicatoron=False, variable=rbRefgenome, value=("hs38d1", 2580))
        rbhs38d1.grid(column=0, row=rowcnt, padx=5, pady=3)

        rbhs38d1s = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38d1s'], font=font['14'],
                                indicatoron=False, variable=rbRefgenome, value=("hs38d1s", 2581))
        rbhs38d1s.grid(column=1, row=rowcnt, padx=5, pady=3)
        rowcnt += 1

        if not inBAM:       # Only display if in BAM (cannot select Unknown if in Align command)
            unknown = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['Unknown'], font=font['14'],
                                  indicatoron=False, variable=rbRefgenome,
                                  value=("Unknown", wgse.BAM.SNCount if inBAM else 0))
            unknown.grid(column=0, row=rowcnt, columnspan=2, padx=5, pady=3)
            rowcnt += 1

        Button(askRefgenWindow, text=wgse.lang.i18n['Done'], font=font['14'],
               command=self.ask_reference_genome_done
               ).grid(column=0, row=rowcnt, columnspan=2, padx=5, pady=3)
        rowcnt += 1

        # This only applies when asking the user for the refgenome (closest) when we could not figure it out
        if inBAM:       # Based on what we already know about BAM, disable some obvious choices that it cannot be
            if not (wgse.BAM.Build == 19 and wgse.BAM.SNTypeC == "Chr"):
                rbhg19.configure(state='disabled')

            if not (wgse.BAM.Build == 37 and wgse.BAM.SNTypeC == "Num"):
                rbhs37d5.configure(state='disabled')
                rbhs37m.configure(state='disabled')
                rbgrch37.configure(state='disabled')

            if not (wgse.BAM.Build == 37 and wgse.BAM.SNTypeC == "Chr"):
                rbhg37.configure(state='disabled')
                rbhg37wg.configure(state='disabled')

            if not (wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Num"):
                rbgrch38.configure(state='disabled')
                rbhs38d1.configure(state='disabled')
                rbhs38d1s.configure(state='disabled')

            if not (wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Chr"):
                rbhs38.configure(state='disabled')
                rbhs38dh.configure(state='disabled')
                rbhg38.configure(state='disabled')

            if not (wgse.BAM.SNTypeC == "Num"):
                rbthgv27.configure(state="disabled")
                rbthgv20.configure(state="disabled")
                rbhppv11.configure(state="disabled")
                rbhppv1.configure(state="disabled")
        else:
            rbgrch37.configure(state='disabled')    # BWA Index is not working on EBI GRCh models so disable
            rbgrch38.configure(state='disabled')

        self.askRefgenWindow = askRefgenWindow
        askRefgenWindow.update()
        askRefgenWindow.grab_set()  # ala mainloop for Toplevel, transient windows
        askRefgenWindow.wait_window()

        # If hit "x' to close window or Done without selecting, then string variable is empty
        refgenome, snCnt = rbRefgenome.get().split() if rbRefgenome.get() else ("Unknown", 0)
        if inBAM:     # Store result in BAM class settings
            wgse.BAM.Refgenome = refgenome
            # wgse.BAM.SNCount = snCnt  # BAM file has its actual SN Count stored which may be different
            if wgse.BAM.Refgenome == "Unknown":
                wgse.BAM.Refgenome = None   # We tried everything; leave it unknown (None)
            if not wgse.BAM.RefMito:
                wgse.BAM.RefMito = "Yoruba" if wgse.BAM.Refgenome == "hg19" else \
                                   "rCRS" if wgse.BAM.Refgenome else \
                                   None
            wgse.BAM.RefgenomeNew = None  # Todo How to set up new refgenome after the fact?
        return refgenome, snCnt

    def ask_reference_genome_done(self):
        try:
            self.askRefgenWindow.destroy()
        except:
            pass


# The rest of this file defines a stand-alone application to process the Final Assembly Reference Genome files
#
#  A Final Assembly (FA) Human Reference Genome Processor to create keys to match to SAM/BAM/CRAM file header
#  Used to create table to figure out what FA file was used to create a BAM; that can then be used for CRAM, etc
#  FA Files are compressed in BGZF format. As we will only read, we can use standard ZIP readers.
#  Historically, .fa/.fna/.fasta files have .gz suffix if compressed; none if plain text
#  Sometimes, marked .gz or otherwise compressed files use .zip and not .bgz (aliased as .gz)
#
#  UPDATE: Have since learned this same information is in the DICT file which is needed to process the FA here.
#          So much of the code can be simplified to simply call samtools dict if that file is missing

#
#  Remaining functionality moved to process_refgenomes.sh for ease of use. May bring back here later.
#
'''
import locale
import gzip
import hashlib
import os.path

from Bio import SeqIO       # Must install Biopython from PIP first

# Table definitions; initialized at end of file.
refgen_SeedSource = {}
refgen_new = {}
reference_genomes = {}

def calc_md5_large_file(oFN, block_size=256*128):
    md5 = hashlib.md5()
    with open(oFN, 'rb') as f:
        for chunk in iter(lambda: f.read(block_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def characterize_ref_genome(bamsq_header):
    model = {'Type': '', 'Who': '', 'Gen': 0, 'Mito': ''}

    # Generally using length of Chr1 (19/37 vs 38) and Chromosome name (1 vs Chr1) to determine reference genome
    # Using purely the length to determine Yoruba vs rCRS (not sure can check for RSRS)
    # Special exception for hs37d5 for Marko :)
    GRCh = bamsq_header.get('1', 0)
    hg   = bamsq_header.get('CHR1', 0)
    Mito = bamsq_header.get('M', bamsq_header.get('MT', bamsq_header.get('CHRM', bamsq_header.get('CHRMT', 0))))

    # Todo handle more than 38, 37, 19 (add 18, 36, etc)
    model['Who']   = 'GRCh' if GRCh else ('hg' if hg else 'Unk') # No Chromosome 1 in model?
    model['Gen']   = 38 if 248956422 in [GRCh, hg] else (37 if 249250621 == GRCh else (19 if 249250621 == hg else 0))
    model['Mito']  = 'rCRS' if Mito == 16569 else ('Yoruba' if Mito == 16571 else 'Unk')
    model['Type']  = 'hs37d5' if bamsq_header.get('HS37D5', 0) > 0 else model['Who'] + str(model['Gen'])
    # DEBUG(f'Model: {model}')
    return model


# DEPRECATED; available from .dict file which is needed anyway. So create that and simply cut from there
# Use local refgenome FA to create a pseudo BAM SN header (upcase, sort -d) to generate a unique md5 hash key
def proc_ref_genome(oFN):
    DEBUG(f"Starting Final Assembly: {oFN}")
    # Process FASTA File; simple 2 line inner loop due to BioPython FASTA file processor
    bamsq_header = {}
    with gzip.open(oFN, "rt") as f:      # BGZF files can be simply unzip'ed by gzip
        for record in SeqIO.parse(f, "fasta"):
            bamsq_header[record.id.upper()] = len(record.seq)

    # Determine characteristics of reference genome model (base model, etc)
    model = characterize_ref_genome(bamsq_header)

    # Create and write pseudo BAM header (upcased, ASCII sorted)
    # Opening file as binary and using encode() prevents the Win10/DOS changing of \n to \r\n
    with open(oFN.replace('.gz','.txt.'), "wb") as outf:
        locale.setlocale(locale.LC_ALL, "C")  # Setting Locale to do ASCII sort always; need "sort -d" to match as well
        pseudo_header = ""
        for key in sorted(bamsq_header, key=locale.strxfrm):
            line = f'SN:{str(key)}\tLN:{bamsq_header[key]}\n'
            outf.write(line.encode()) # Keep in UNIX format; write out for human documentation / validation
            pseudo_header += line
        locale.setlocale(locale.LC_ALL, '')   # Restore back to local, default, for print formatting below

    # Get md5 hash of pseudo BAM header
    md5hash = hashlib.md5(pseudo_header.encode())

    # Create table entry for Reference_Genome Known Models
    model_WhoGen = model['Who'] + str(model['Gen'])
    BN = os.path.basename(oFN)
    BN_size_KB   = round(os.path.getsize(oFN)) / 1000    # In KB, purposely, to avoid OS file block dependencies (not 1024, matching to Unix ls -l command)
    BN_bamhd_md5 = md5hash.hexdigest()
    BN_file_md5  = calc_md5_large_file(oFN)
    print(f"'{BN_bamhd_md5}': ['{model['Type']}',\t'{model_WhoGen}',\t'{model['Mito']}',\t'{BN}',\t'{BN_size_KB:n} KB',\t'',\t'{BN_file_md5}'],")
    desc, URL1, URL2, URL3 = refgen_SeedSource.get(BN, ['', '', '', ''])
    reference_genomes[BN_bamhd_md5] = [model['Type'], model_WhoGen, model['Mito'], BN, BN_size_KB, BN_file_md5, desc, URL1, URL2]


if __name__ == '__main__':
    wgse.DEBUG_MODE = True
    RefGen_oFP = "D:\\Reference_Genomes\\"

    # Bio.seqIO solution to read FA file
    for file in refgen_SeedSource:
    #    put up Please Wait (30 minutes)
        proc_ref_genome(RefGen_oFP + file)
    #    cancel wait()
    print(reference_genomes)
'''
''' # PyFaidx solution -- not any faster than using Bio.SeqIO
    for file in dir_list_of_files:
        if file has .fai or file is bgzf compressed:
            put up Please Wait (30 minutes) if no .fai as must create
            process with pyfaidx
            cancel wait() if started
        else:   # My code is just as efficient using SeqIO from BIO
            proc_ref_genome(path+file)
    with open(RefGen_oFP + "Reference_Genomes.txt", 'w') as f:
        f.write(str(reference_genomes))
'''
'''
    # import pybam    # Written in Python2; having trouble converting
    # for bam in pybam.read("I:\WGS\RandyYseq\\22731_bwa-mem_hg19_sorted.bam"):
    #    print(bam.file_chromosome_lengths)  # quick BAM header without running bash script and samtools view

    # PySAM solution -- part of Bioconda; uses cython and htslib library to give Samtools.BCFtools direct access
    # But install requires C compiler; so not likely sustainable on Win10 where Bioconda not available.
'''
################################################################################################################
# Reference Tables
#
#  Really need to be moved to external files and thus maintained independently of this code. None currently used anyway!
#
# Key table of reference genomes (AUTO GENERATED FROM ABOVE CODE -- DO NOT EDIT DIRECTLY)
#  BAMHeader MD5 Sum: Jon Rhys MD5Sum of the SQ records in the header; used as key to entry here. Taken from SN and Len
#  Generic Model: HG19/38 has chr1-chr22 names, GRCh37/38 has 1-22 names
#  Mito Model: Which mitochondrial reference model is used (should be able to detect but just for documentation here)
#  File Name: The local file name of the model to supply to the tools (FBS form)
#  File Size:
#  URL:  The URL to do a wget/curl to get the file
#  RefGenFile MD5 Sum: MD5Sum of actual reference genome file name (to validate download copy)
reference_genomes = {
'ee4efe40ebd6f9468dab89963dcc5b65': ['hg19',	'hg19',	'Yoruba',	'hg19.fa.gz',	'948,731 KB', '',	'806c02398f5ac5da8ffd6da2d1d5d1a9'],
'eed540239f34de0b80734add1779cb5b': ['hg19',	'hg19',	'Yoruba',	'hg19.p13.plusMT.fa.gz',	'979,165 KB', '',	'7707462fc100c7d987c075bc146b16ae'],
'8ed6919c076595f09e1028325dc86fbe': ['hg19',	'hg19',	'rCRS',	'hg19.p13.plusMT.full_analysis_set.fa.gz',	'900,428 KB', '',	'68ba6e8b063e46bbf2eff013d8becadf'],
'5a923ccfc7ee30aebe864817653d41d8': ['hg19',	'hg19',	'rCRS',	'hg19.p13.plusMT.no_alt_analysis_set.fa.gz',	'863,259 KB', '',	'd8cae8ef0f24723a31f224ec9a74f874'],
'9b3bf038da5fef087723525d5174f0bd': ['hg19',	'hg19',	'rCRS',	'hg19_yseq.fa.gz',	'970,288 KB', '',	'09570164f73fbe82a30428fbfdeb4e44'],
'50093bf43ed53848c9c7696752e65c2f': ['hg19',	'hg19',	'rCRS',	'hg19_wgse.fa.gz',	'881,223 KB', '',	'4edc36f1c0c7b854d48853197cb7dd1e'],
'bd894134bddba260df88a90123a2ee9c': ['hg38',	'hg38',	'rCRS',	'hg38.fa.gz',	'983,659 KB', '',	'1c9dcaddfa41027f17cd8f7a82c7293b'],
'f5abbfc73ccdb2e8356003d230bce0ce': ['hg38',	'hg38',	'rCRS',	'hg38p12.fa.gz',	'998,672 KB', '',	'242559a5aac0e1720a6c3396674f36ad'],
#                                                                'human_g1k_v37.fasta.gz'
'88b5f4bbaaec4a097b3aa160bf776093': ['GRCh37',	'GRCh37',	'rCRS',	'human_g1k_v37_WGSE.fasta.gz',	'882,999 KB', '',	'd3328225a911adcf3593e81721793819'],
#                                                                'hs37d5.fa.gz'
'92c90b84dac3d9d8402bf2f2dbe95961': ['hs37d5',	'GRCh37',	'rCRS',	'hs37d5.fa.gz',	'892,326 KB', '',	'd10eebe06c0dbbcb04253e3294d63efc'],
'f8812b68f2aa4aa592e742e081f82807': ['GRCh37',	'GRCh37',	'rCRS',	'Homo_sapiens.GRCh37.latest.dna.toplevel.fa.gz',	'993,582 KB', '',	'a4803af89b6c7a9a0db87aa9f8e73461'],
'abf297943fca66859d781349e8e28687': ['GRCh38',	'GRCh38',	'rCRS',	'Homo_sapiens.GRCh38.99.dna.toplevel.fa.gz',	'1.10765e+06 KB', '',	'ae85c8481141a9399ec7576ad6b0e259'],
'5c38f3f6774ece1a75c4d5d55010ced7': ['hg19',	'hg19',	'rCRS',	'GRCh37.p13.genome.fa.gz',	'804,606 KB', '',	'490904ef2561472d0801c704ffc85761'],
'a08daf6f9f22170759705fd99e471b62': ['hg38',	'hg38',	'rCRS',	'GRCh38_full_analysis_set_plus_decoy_hla.fa.gz',	'918,931 KB', '',	'9513ce08c458ac88f8411dcf01097a1f'],
'42c06a6446741565da8cae7bbeb67184': ['hg38',	'hg38',	'rCRS',	'GRCh38.p13.genome.fa.gz',	'889,409 KB', '',	'1b4b49100b0df70e097a4ab9f4db8b70'],
'236b4fa639deaaa51691eb16ea92a0f9': ['hg19',	'hg19',	'rCRS',	'GCA_000001405.14_GRCh37.p13_full_analysis_set.fna.gz',	'900,426 KB', '',	'9a705a08e567abbcfd6bfc7a0c6d87cd'],
'0b94854439441bc45d1b061eedb2ac0e': ['hg19',	'hg19',	'rCRS',	'GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz',	'863,259 KB', '',	'21e7083e341404f752682e3b289ba454'],
'53f744f1771eabfded06778564b931bc': ['hg38',	'hg38',	'rCRS',	'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz',	'903,027 KB', '',	'92d3088bc1c97478c37c2e13b620d50c'],
#                                                               'GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz'
'c9ef5f042a951d511567dd4ea9bd126e': ['hg38',	'hg38',	'rCRS',	'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz',	'872,950 KB', '',	'a08035b6a6e31780e96a34008ff21bd6'],
'ee77d771adff64d7ee6357316ee24b6b': ['hg38',	'hg38',	'rCRS',	'GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz',	'874,638 KB', '',	'a056c57649f3c9964c68aead3849bbf8']
}


# New approach.  Use SN count instead of MD5 on SN and LN fields in BAM header. Seems possibly more robust / reliable.
# This after discovering the simple approach of model identification via Chr1 LN and SN name (with/out chr prepend) is
# not accurate to determine the model type.  GCA_000*fna.gz has "chrNN" but is labeled in v2 as GRCh38.
# In python if multiple dictionary entries with the same key are specified, only the last is saved.  We have avoided
# doing a double level dictionary for the few duplicate keys there are.  May be an issue later and cause is to
# readdress using the MD5 solution key solution.
refgen_new = {
     25: [True,  "hg19_WGSEv2.fa.gz", "chrM"],  # WGSE v2 hg19  (really in error, is primary only)
     84: [True,  "human_g1k_v37.fasta.gz", "MT"],  # WGSE v2 GRCh37 (base 1K Genome analysis model)
     85: [False, "GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz", "chrM"],  # aka hs38.fa.gz, build 37 patch 13
     86: [True,  "hs37d5.fa.gz", "MT"],  # WGSE v2 hs37d5
     93: [False, "hg19.fa.gz", "chrM"],  # TRUE hg19 file (not Marko's WGSEv2 short one), no patches, Yoruba (only one)
    195: [True,  "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz", "chrM"],  # WGSE v2 GRCh38 ; aka hs38.fa.gz
    297: [False, "Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz", "MT"], # EBI Ensemble release 75 (Build 37, patch 13)
#   297: [False, "GRCh37.p13.genome.fa.gz", "chrM"],
    298: [False, "GCA_000001405.14_GRCh37.p13_full_analysis_set.fna.gz", "chrM"], # 1K Genome Build 37 patch 13
    389: [False, "Homo_sapiens.GRCh38.89.dna.toplevel.fa.gz", "MT"],  # EBI/Ensemble release 89 (patch 2?)
    455: [False, "hg38.fa.gz", "chrM"],  # Initial model, no patches
    456: [False, "GCA_000001405.15_GRCh38_full_analysis_set.fna.gz", "chrM"],  # Companion to 195
    524: [False, "Homo_sapiens.GRCh38.86.dna.toplevel.fa.gz", "MT"],  # EBI/Ensemble 86/87 (original and patch 1)
    555: [False, "Homo_sapiens.GRCh38.91i.dna.toplevel.fa.gz", "MT"],  # EBI/Ensemble 91i; 90i, 88 also same count
    593: [False, "Homo_sapiens.GRCh38.97i.dna.toplevel.fa.gz", "MT"],  # EBI/Ensemble 97i; 96i-92i also same count
    595: [True,  "hg38p12.fa.gz", "chrM"],  # WGSE v2 hg38  (Full Analysis Set, patch 12)
    639: [False, "Homo_sapiens.GRCh38.current_fasta.dna.toplevel.fa.gz", "MT"], # EBI/Ensemble aka 100i; 99i and 98i
#   639: [False, "GRCh38.p13.genome.fa.gz", "CHRM"], # Note: unique upcase CHRM nomenclature
   2580: [False, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz", "chrM"],  # hs38DH no alt
   2841: [False, "GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz", "chrM"],  # hs38DH Full
   3366: [False, "GRCh38_full_analysis_set_plus_decoy_hla.fa.gz", "chrM"]  # Full Analysis Set plus decoy HLA
}


# This is a HAND-CRAFTED seed table of URLs for specific local files. Local file name is dictionary key.
# Note: the original file name may conflict (non-unique) and have to be renamed.  Also may not be compressed or .zip compressed.
refgen_SeedSource = {
#   Key (Local FN)   Description           Primary URL                 Secondary URL                   Tertiary URL
    'hg19.fa.gz':
        ['UCSC HGP HG19 Original Release (Nov 2009)',
         'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz', '', ''],
    'hg19.p13.plusMT.fa.gz':
        ['UCSC HGP HG19 Patch 13 (Jan 2020) (latest)',
         'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz', '', ''],
    'hg19.p13.plusMT.full_analysis_set.fa.gz':
        ['UCSC HGP HG19 Path13 Full Analysis)',
         'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.full_analysis_set.fa.gz',
         '', ''],
    'hg19.p13.plusMT.no_alt_analysis_set.fa.gz':
        ['UCSC HGP HG19 Patch 13 No Alt Analysis',
         'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz', '', ''],
    'hg19_yseq.fa.gz':
        ['UCSC HGP HG19 version? (ySeq.net disto)',
         '', '', ''],
    'hg19_wgse.fa.gz':
        ['UCSC HGP HG19 version? (WGSE v2 disto)',
         'WGSE', '', ''],
    'hg38.fa.gz':
       ['UCSC HGP HG38 Original Release',
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz', '', ''],
    'hg38p12.fa.gz':
        ['UCSC HGP HG38 Path 12 (Aug 2018) (latest)',
         'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz', '', ''],
    'hg38_yseq.fa.gz':
        ['UCSC HGP HG38 version? (ySeq.net disto)',
         '', '', ''],
    'hg38_wgse.fa.gz':
        ['UCSC HGP HG38 version? (WGSE v2 disto)',
         'WGSE', '', ''],
#    'human_g1k_v37.fasta.gz':              # line 411, gzip.py, OSError: Not a gzipped file (b'\x01w') (deep)
#        ['1K Genome Project Phase 1 GRCh37',
#         'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz', '', ''],
    'human_g1k_v37_WGSE.fasta.gz':
        ['1KGenome Project Phase 1 GRCh37 (WGSE v2 disto)',
         'WGSE', '', ''],
#    'hs37d5.fa.gz':                        # line 411, gzip.py, OSError: Not a gzipped file (b'\x01') (deep)
#        ['1K Genomes hs37d5 (ncbi)',
#         'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
#         '', ''],
    'hs37d5_unknown_source.fa.gz':
        ['1K Genomes hs37d5 (unknown source)',
         '', '', ''],
    'hs37d5_WGSE.fa.gz':
        ['1K Genomes hs37d5 (WGSE v2 disto)',
         'WGSE', '', ''],
    'Homo_sapiens.GRCh37.latest.dna.toplevel.fa.gz':
        ['EBI/Ensemble GRCh37 (latest) (Genbank GCA_1405.14)',
         'ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz', '', ''],
    'Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz':
        ['EBI/Ensemble GRCh37 p13 (release 75) (Genbank GCA_1405.14)',
         'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz', '', ''],
    'Homo_sapiens.GRCh38.99.dna.toplevel.fa.gz':
        ['EBI/Ensemble GRCh38 p13 (release 99) (Genbank GCA_1405.28)',
         'ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz', '', ''],
    'GRCh37.p13.genome.fa.gz':
        ['',
         'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz', '', ''],
    'GRCh38_full_analysis_set_plus_decoy_hla.fa.gz':
        ['1KGenome Project Phase 3 GRCh38 Full Analysis plus decoy hla (Jul2015)',
         'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa',
         '', ''],
    'GRCh38.p13.genome.fa.gz':
        ['EBI/Ensemble 1KGenome GRCh38 patch 13 (Gencode release 33, Ensemble 99) (Jan2020)',
         'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.p13.genome.fa.gz', '', ''],
    'GRCh38.p13.Rel32.fa.gz':
        ['EBI/Ensemble 1KGenome GRCh38 patch 13 (Gencode release 32, Ensemble 99)',
         'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz', '', ''],
    'GCA_000001405.14_GRCh37.p13_full_analysis_set.fna.gz':
        ['NCBI GenCode GRCh37 Patch 13 Full Analysis (GCA_1405.14)',
         'ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/seqs_for_alignment_pipelines/GCA_000001405.14_GRCh37.p13_full_analysis_set.fna.gz',
         '', ''],
    'GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz':
        ['NCBI GenCode GRCh37 Patch 13 No Alt Analysis (GCA_1405.14)',
         'ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/seqs_for_alignment_pipelines/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz',
         '', ''],
    'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz':
        ['NCBI GenCode GRCh38 Full Analysis (GCA_1405.15)',
         'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz',
         '', ''],
#    'GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz':    # File "gzip.py", line 482, EOFError: Compressed file ended before the end-of-stream marker was reached
#        ['NCBI GenCode GRCh38 Full Analysis plus HS38D1 (GCA_1405.15)',
#         'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz',
#         '', ''],
    'GCA_000001405.15_GRCh38_no_alt_analysis_set_WGSE.fna.gz':
        ['NCBI GenCode GRCh38 No Alt (WGSE v2 disto)',
         'WGSE', '', ''],
    'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz':
        ['NCBI GenCode GRCh38 No Alt Analysis (GCA_1405.15)',
         'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz',
         '', '']
#    'GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz':  # File "C:\Program Files\JetBrains\PyCharm Community Edition 2019.3.4\plugins\python-ce\helpers\pydev\pydevd.py", line 1434, in _exec
#        ['NCBI GenCode GRCh38 No Alt Analysis plus HS38D1 (GCA_1405.15)',
#         'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna,gz',
#         '', '']
}


