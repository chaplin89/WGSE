# coding: utf8
#
# Reference Library Class subsystem
#
# Part of the
# WGS Extract (https://wgse.bio/) subsystem
#
# Copyright (C) 2020-2023 Randy Harr
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
import csv
import os.path
import shutil

import settings as wgse
from commandprocessor import run_bash_script
from utilities import DEBUG, nativeOS, is_legal_path, unquote, wgse_message


class ReferenceLibrary:
    """
    In preparation for a full function Reference Library capability (library maintenance of multiple types of objects;
    manipulation, automatic keep-up-to-date, let users add entries, etc), creating this class now to start.
    For now, mostly to capture and handle the Reference Library file-name / access variables that were in settings.
    Too much functionality was buried in mainwindow that needed to be pulled out to handle the settings storage and
    restore.  So better to put it into a class here than mainwindow which should be UI focused.
    """

    # Todo additional Liftover files: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ and
    #  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ should be added to reference library and defined
    #  for use as should VCF annotations, BED files, etc

    def __init__(self, default_FP):
        self.valid = False      # True if current reference library directory exists and has content
        self.set = None         # True if reflib dir is not default; has been overriden
        self.default_FP = None  # Default Reference Library location FP (initial call; back to if unset)
        self.FP  = None         # Reference Library File (Directory) Pointer
        self.oFP = None
        self.gen_FP  = None     # (Reference) Genome File (Directory) Pointer
        self.gen_oFP = None
        self.FB  = None         # Final (leaf) directory name of Reference Library location (e.g. "reference")
        self.cma_FP = None      # Chip MicroArray Templates and support File (Directory) Pointer
        self.cma_oFP = None
        self.liftover_chain_oFN = None      # Liftover Chain pointer (currently Build38 -> Build37 only)
        self.askRefgenWindow = None

        # Index to list of each row in genomes.csv file in new dictionary; first two entries of file not there
        # IMPORTANT: this order defines columns in genomes.csv file and must match it (no column labels in file)
        # self._sortindex = -2      # Drop first colunn; used only by user editing the genomes.csv file
        # self._refgenome = -1      # Second column becomes the dictionary key with remaining columns a list entry
        self._sn_count = 0
        self._sn_naming = 1
        self._build = 2
        self._class = 3
        self._source = 4
        self._final_file = 5
        self._init_file = 6
        self._url = 7
        self._menu_cmd = 8
        self._descr = 9
        self._md5f = 10
        self._mdfi = 11
        self._sizef = 12
        self._sizei = 13

        # Genomes.csv column definitions in post dictionary creation form (per genome)
        self._keys = []             # Keys used in dictionaries (complete ordered set from genomes.csv) for convenience
        self._genomes = {}          # genomes.csv in dictionary form using _refgenome as key into dictionary

        # Lists of genomes to display to the user (more desciptive; composed of multiple _genomes() entries in a string
        self._askgenomes = {}       # List of genome descriptions to display to user when selecting a ref genome
        self._availist = []         # List of ref genome descriptionss available at that time for user to select
        self._revkey = []           # Key into _genomes for corresponding _availist list entry

        self.change(default_FP)     # Call Reference Library location Change with Default location (for initial setup)

    def change(self, dir_FP):
        """ Change ther reference library folder location. Used to set initial one if not restored from settings. """

        # Directory must not be the root file system (also check for null value early on)
        if not dir_FP or not is_legal_path(dir_FP):
            if wgse.lang.i18n:
                wgse_message("error", 'InvalidRefLibTitle', True,
                             wgse.lang.i18n['errRefLibPath'].replace("{{DIR}}", dir_FP))
            else:
                DEBUG(f'***ERROR: Reference Library directory path is not valid: "{dir_FP}"')
            return

        dir_FP += '/' if dir_FP[-1] != '/' else ''          # Assure trailing slash; universalOS so forward slash
        if os.path.basename(dir_FP[:-1]) != "reference":    # Always name the directory "reference/" (end of path)
            dir_FP += "reference/"

        # No change. Silently return
        if dir_FP == self.FP:
            return

        dir_oFP = nativeOS(dir_FP)

        gen_FP  = f'{dir_FP}genomes/'
        gen_oFP = nativeOS(gen_FP)
        cma_FP = f'{dir_FP}microarray/'
        cma_oFP = nativeOS(cma_FP)

        # Does new reflib directory exist and have content ?
        valid = os.path.isdir(dir_oFP) and os.path.isdir(gen_oFP) and os.path.isdir(cma_oFP)

        # Initial setup of Default dir in WGS installation folder (reference/); may not exist if moved in settings
        if not self.default_FP and not valid:
            DEBUG(f'***WARNING: the default Reference Library does not exist or is empty (moved?): "{dir_FP}"')

        elif self.default_FP and not valid and not self.valid:
            if wgse.lang.i18n:
                wgse_message("error", 'InvalidRefLibTitle', True,
                             wgse.lang.i18n['errRefLibEmpty'].replace("{{DIR}}", dir_oFP))
            else:
                DEBUG(f'***WARNING: Existing and new Reference Library have no content: "{self.FP}" "{dir_FP}"')
            return

        # Does current reflib setting exist and have content but new one does not? If so, move current content to new
        if not valid and self.valid:

            # Get the size of the existing reference library; format into a usable string for the user
            sz = sum(os.path.getsize(os.path.join(dirpath, filename))
                     for dirpath, dirnames, filenames in os.walk( self.oFP) for filename in filenames)
            scaled, mod = (sz/10**6, "MB") if sz > 10**6 else (sz/10**3, "KB") if sz > 10**3 else (sz, "")
            sz_str = f'{scaled:,.0f} {mod}'

            # Query the user to go ahead with the move; if UI and lang setup (must hit OK button)
            if wgse.lang.i18n:
                mesg = wgse.lang.i18n['ConfirmationRefLib'].replace(
                    "{{SIZE}}", sz_str).replace("{{SRC}}", self.FP).replace("{{DST}}", dir_FP)
                if not wgse_message("okcancel", 'ConfirmationRefLibTitle', True, mesg):
                    return
            DEBUG(f'*** WARNING: Moving {sz_str} from folder {self.FP} to {dir_FP}.')

            # Perform the move of the existing reference library to the new location
            shutil.move(self.oFP, os.path.dirname(dir_oFP[:-1]), copy_function=shutil.copy2)

        # Now that we have error checked (confirmed) values, go ahead and replace values in the Class instance
        if not self.default_FP:                     # If first time initializsing, then set the default
            self.default_FP = dir_FP
        self.set = dir_FP != self.default_FP        # Meaning set by user or settings; not default
        self.FP = dir_FP        ;  self.oFP = dir_oFP
        self.gen_FP = gen_FP    ;  self.gen_oFP = gen_oFP
        self.cma_FP = cma_FP    ;  self.cma_oFP = cma_oFP
        self.FB = os.path.basename(dir_FP[:-1])     # Should always be "reference"
        self.valid = os.path.isdir(dir_oFP) and os.path.isdir(gen_oFP) and os.path.isdir(cma_oFP)

        self.liftover_chain_oFN = f'{self.oFP}hg38ToHg19.over.chain.gz'
        # Todo change from static value to method call so can have chain files for multiple liftover directions

        self._load_genomes_csv()        # May have changed source of genomes.csv file

        DEBUG(f'New Reflib Directory: {self.FP}')

    @staticmethod
    def getmt(oFN):
        return os.path.getmtime(oFN) if os.path.exists(oFN) and os.path.isfile(oFN) else 0

    def _load_genomes_csv(self):
        """ Initialize the genomes directory from the genomes.csv file. Create a dictionary of lists.
            NOTE: genomes.csv has two entries for the PythCode column (i.e. the ref_genome internal unique keyword).
            The two entries should be identical except for the source name, download URL, Library Menu Label,
            and (longer) Description. We do not really care about those 2nd entries as we do not use those alt entries
            in python code.  We call get_refgenomes.sh that reads in genomes.csv itself. Otherwise, we
            would need to use the sort # field as the unique dictionary key. And then go through the array for
            each lookup.  For now, we simply overwrite the first entry with the second. The PythCode is then
            unique to each ref_genome and used as the dictionary key.
            FURTHER NOTE: genomes.csv is missing most of the patched models (p1-13/14 of EBI 37 / 38). So the
            reference model for BAMs dependent on those models will still not be found.
        """
        # Let's first make sure we have the latest and best genomes.csv file (see zcommon.sh genomes.csv loader also)
        # Todo check against latent online and use if newer than programseed
        # Todo because check everytime, do not copy; just use current best in situ; in case user changes local copy
        programseed = wgse.prog_oFP + 'seed_genomes.csv'        # Now delivered with program package in program/ folder
        reflibseed = self.oFP + 'seed_genomes.csv'              # Initially delivered in reflib package in reference/
        genomescsv = self.gen_oFP + 'genomes.csv'               # Only used in reference/genomes/ folder

        programseed_mt = self.getmt(programseed)                # File modified time in program/ folder (if exists)
        reflibseed_mt = self.getmt(reflibseed)                  # File modified time in reference/ folder (if exists)
        genomescsv_mt = self.getmt(genomescsv)                  # File modified time in reference/genomes folder

        if reflibseed_mt > genomescsv_mt:
            shutil.copy2(reflibseed, genomescsv)        # Copy while retaining attributes such as modified time
            genomescsv_mt = self.getmt(genomescsv)      # Update modified time for reference/genomes/genoms.csv file

        if programseed_mt > genomescsv_mt:
            shutil.copy2(programseed, genomescsv)
            genomescsv_mt = self.getmt(genomescsv)           # Update modified time for reference/genomes/genoms.csv file

        if genomescsv_mt == 0:
            DEBUG(f'*** ERROR No valid genomes.csv file found at {genomescsv}')
            return 1

        # With latest, valid genomes.csv file, read in content
        with open(genomescsv) as genomesf:
            genomes = csv.reader(genomesf)
            next(genomes)       # Skip first title / comment row

            self._keys = []
            for genome in genomes:
                # DEBUG(f'{genome}')

                # Convert "." shortcut in genomes.csv final file name column to have initial file name value
                # note: working from genomes.csv entry directly before converting to dictionary entry; so +2
                if genome[self._final_file+2] == '.':
                    genome[self._final_file+2] = genome[self._init_file+2]

                key = genome[1]     # column 1 is unique and the key into the dictionary
                cur = genome[2:]    # Skip column 0 (sort order) and column 1 (key); load the rest as list

                # Implementing an ordered set using a list to mimic an ordered dictionary
                if self._genomes.get(key, "missing") != "missing":      # Remove existing key before adding again
                    self._keys.remove(key)
                self._keys += [key]         # This is the list of refgenome codes used inside the program
                self._genomes[key] = cur    # Note: overwrites same key entries (NIH vs EBI) so ignore init elements

                # Construct the display string to the user for any drop-down menu to select a reference genome
                #  key: _sn_count, _sn_naming, _build, _class, _source, _final_file, _init_file,
                #       _url, _library_menu_cmd, _descr, _md5f, _md5i, _sizef, _sizei
                title = cur[self._menu_cmd].replace(" (*)", "").replace(" 3x", "")  # Clear Library menu items in name
                if not any(s in cur[self._class] for s in ["hs37h", "hs37fh"]):     # Todo Why not every entry?
                    title = title.replace(" (NIH)", "").replace(" (EBI)", "")       # Remove model source component

                name = cur[self._sn_naming]
                mito = "Yoruba" if "Yoruba" in title else "rCRS"
                sns = cur[self._sn_count]

                self._askgenomes[key] = f'{title:>30}, {name:>5}, {mito:>6}, {int(sns):5,d} SNs'     # 56 char field

        DEBUG(f'Loaded genomes.csv, {len(self._genomes)} entries')

    def refgenome_codes(self):
        return self._keys

    def available_refgenomes(self, inBAM=True):
        """ Create a combobox (dropdown list) of possible refgenomes to select from the loaded genoms.csv file.
            Combobox only takes a list of strings that it displays. It then returns one of those strings as the
            selected result. So we have to use the constructed dictionary "_askgenomes[]" which has the display
            strings and the refgenome code as its key. Then use the returned string to look up the key in there.
            inBam specifies that the list should be trimmed to only those entries that are valid for the BAM
            settings. We are asking the user to select the closest ref genome to what they see is in the BAM.
            So impossible matches like wrong chromosome names or major build are trimmed out. And an "Unknown"
            added so a user can select that as well. Otherwise, the whole list from genomes.csv is returned.
        """

        self._availist = []
        self._revkey = []
        for key, genome in self._genomes.items():

            # If inBAM, then prune out genomes that do not match the chr naming and major build of the BAM
            if inBAM and wgse.BAM and \
              (genome[self._sn_naming] not in [wgse.BAM.SNTypeC, wgse.BAM.SNTypeM] or
               wgse.BAM.Build != genome[self._build]):
                # or genome[self._sn_count] < wgse.BAM.SNCount
                continue

            # Create list with composed string of each reference genome to use in drop down (and corresponding key)
            self._availist += [self._askgenomes[key]]
            self._revkey += [key]

        # For user choosing closest refgenome for BAM, allow selecting "unknown" (add as last entry in dictionary)
        if inBAM:
            self._availist += [f'{"Unknown":>24}']
            self._revkey += ["unknown"]

        return self._availist

    def key_from_availist(self, refstring):
        """
        When a user selects a refgenome from the drop down list, we are returned with that entry of the availist array.
        So use that string to look into the _availlist to get an index back into the _revkey for the refgenome code.
        """

        # Note: returns "unknown" if so selected; is not a key into _genomes though
        return next(self._revkey[i] for i, askref in enumerate(self._availist) if askref == refstring)

    def get_refgenome_qFN(self, ref_genome):
        """  Simple lookup and return of the (final) file name for the reference genome key """

        if ref_genome and ref_genome != "unknown":
            genome = self._genomes.get(ref_genome, [])
            if len(genome) > 0:
                return f'"{self.gen_FP}{genome[self._final_file]}"'
        # else return null

    def get_refgenome(self, refgenome_qFN):
        """ Simple reverse lookup to return of the refgenome key given the quoted refgenome final file name """
        return next(key for key, entry in self._genomes.items()
                    if refgenome_qFN.endswith(entry[self._final_file] + '"'))

    def installed_refgenome(self, ref_genome="", refgenome_qFN=""):
        installed = True

        if ref_genome and not refgenome_qFN:
            refgenome_qFN = self.get_refgenome_qFN(ref_genome)

        if not refgenome_qFN:
            return not installed

        refgenome_oFN  = nativeOS(unquote(refgenome_qFN))
        wgse_oFN = refgenome_oFN.replace("fa.gz", "wgse").replace("fna.gz", "wgse").replace("fasta.gz", "wgse")

        refgenome_mt = self.getmt(refgenome_oFN)    # Note: returns 0 if file does not exist
        wgse_mt = self.getmt(wgse_oFN)              # Note: returns 0 if file does not exist

        if not ref_genome:
            ref_genome = self.get_refgenome(refgenome_qFN)
        target_size = int(self._genomes[ref_genome][self._sizef])

        # Standard case -- refgenome file found, correct size and processed already
        if 0 < refgenome_mt < wgse_mt and os.path.getsize(refgenome_oFN) == target_size:
            return installed

    @staticmethod
    def determine_refmodel(Header, Build, RefMito, SNCount, SNTypeC, SNTypeM):
        """
        Determine the reference model used to create a SAM / BAM / CRAM.  Based on Major / Minor / Class genome study
        in https://bit.ly/34CO0vj. Used to rely on SNCount.  But oddball, unrecognized ref genomes in BAMs may have
        same SNcount. So now use key items passed in beside the header. Reference genomes returned here must
        match those in the first column of seed_genomes.csv and cover all entries completely to have them recognized.
        Todo finally implement lookup using the MD5Sum of BAM Header SQ fields, Make use of .wgse files created ...
        Todo when compiling a refgenome. Allow new entries by the user to be recognized and processed here and elsewhere
        """

        if Build == 99:  # Used chr1, chrX or chrY lengths to determine Build99 earlier
            # Use X and Y length to determine T2T / HPP model
            if "LN:62456832" in Header and "LN:154343774" in Header:  # CHM13 v1.1 & HG002 v2 XY
                Refgenome = RefgenomeNew = "THGv20"
            elif "LN:62460029" in Header and "LN:154349815" in Header:  # CHM13 v1.1 & HG002 v2.7 XY
                Refgenome = RefgenomeNew = "THGv27"
            elif "LN:62480187" in Header and "LN:154434329" in Header:  # HG01243 v3 PR1 "Puerto Rican"
                Refgenome = RefgenomeNew = "THG1243v3"
            elif "LN:57227415" in Header and "LN:154259566" in Header:  # CHM13 v1.1 & GRCh38 Y
                Refgenome = RefgenomeNew = "HPPv11"
            elif "LN:57227415" in Header and "LN:154259625" in Header:  # CHM13 v1 X & GRCh38 Y
                Refgenome = RefgenomeNew = "HPPv1"
            elif "LN:62456832" in Header and "LN:156040895" in Header:  # GRCh38 w/ HG002 v2 Y  (ySeq)
                Refgenome = RefgenomeNew = "THGySeqp"
            elif "LN:62460029" in Header and "LN:154259566" in Header:  # CHM13 v1.1 & HG002 v2.7 Y
                Refgenome = RefgenomeNew = "T2Tv20"
            elif "LN:154259566" in Header:  # Must be plain CHM13 v1.1 as special Y not found
                Refgenome = RefgenomeNew = "T2Tv11"
            elif "LN:154259625" in Header:  # Must be plain CHM13 v1 as special Y not found
                Refgenome = RefgenomeNew = "T2Tv10"
            elif "LN:154259664" in Header:  # Must be plain CHM13 v0.9
                Refgenome = RefgenomeNew = "T2Tv09"
            else:  # Think it is T2T / HPP build but cannot recognize Y or X length
                DEBUG(f'Unrecognized T2T / HPP model; no known X and/or Y length in BAM.')
                Refgenome = RefgenomeNew = None
            if RefgenomeNew == "T2Tv20" and SNTypeC == "Acc":
                Refgenome = RefgenomeNew = "T2Tv20a"

        elif SNCount in [85, 86] and Build == 37 and "SN:NC_007605" in Header:
            # 1k/hs37 class have SN:NC007605 (EBV in Numeric naming) sans human_g1k / GRCh37.primary_assembly models
            Refgenome = "hs37d5" if "@SQ\tSN:hs37d5" in Header else "hs37g"
            RefgenomeNew = "1k37" + "g" if SNTypeC == "Num" else "h"

        elif SNCount in [85, 298] and Build == 37 and "SN:chrEBV" in Header and SNTypeC == "Chr":
            # hs37 class with Chr naming and no / full alt analysis sets (1K models @ Genbank; ChrMT odd at UCSC)
            Refgenome = "hs37h" if SNCount == 85 else "hs37fh"
            RefgenomeNew = "1k37h"  # g form does not exist
            if SNTypeM == "chrMT":
                # Special exception for UCSC oddball model hg19.p13.plusMT.*_analysis_set.fa.gz
                Refgenome += "t"
                RefgenomeNew += "t"

        elif SNCount == 84:  # human_g1k (if Num), GRCh37.primary_assembly.genome (if Chr)
            Refgenome, RefgenomeNew = \
                ("hs37-", "1K37g")    if SNTypeC == "Num" else \
                ("GRCh37-", "EBI37h") if SNTypeC == "Chr" else \
                (None, None)

        elif SNCount in [195, 456, 2580, 2581, 2841, 3366] and Build == 38 and \
            ("SN:chrEBV" in Header or "SN:EBV" in Header):
            # hs38DH is 3366 SN count and has SN:HLA- unique; hs38 is 195 and hs38a is 456; all uniquely have chrEBV
            # hs38d1s is sequencing.com model made by hs38d1 with 22_KI270879v1_alt added
            Refgenome, RefgenomeNew = \
                ("hs38DH", "1k38h")   if SNCount == 3366 and "SN:chr22_KI270879v1_alt" in Header else \
                ("hs38d1a", "1k38h")  if SNCount == 2841 else \
                ("hs38d1s", "1k38pg") if SNCount == 2581 and "SN:22_KI270879v1_alt" in Header else \
                ("hs38d1", "1k38h")   if SNCount == 2580 else \
                ("hs38a", "1k38h")    if SNCount == 456 else \
                ("hs38", "1k39h")     if SNCount == 195 else \
                (None, None)

            if SNCount == 2580 and SNTypeC == "Chr" and \
                ("M5:a491618313b78cdca84ae9513e4f4844" in Header or "Verily" in Header):
                Refgenome = "hs38d1v"    # Special Google Verily model; cannot be distinguished without M5 field;
                RefgenomeNew = "1k38ph"  # has very different m5 signature on each SN than hs38d1

        elif SNCount in [93, 297, 455, 639]:  # SNCounts of 6 hg/ebi models; 2 duplicated
            # Handle the hg (UCSC) and GRCh (EBI) models here
            Refgenome = "hg" if SNTypeC == "Chr" else "GRCh"  # 297 & 639 are Patch 13 GRCh models
            Refgenome += str(Build)  # Build already takes into account mito for 19/37 differentiation
            RefgenomeNew = Refgenome.replace("GRCh", "EBI") + ("h" if SNTypeC == "Chr" else "g")

        elif SNCount == 194 and Build == 38 and SNTypeC == "Chr":
            Refgenome = "GRCh38-"
            RefgenomeNew = "hg38h"

        elif SNCount == 25 and Build == 37 and RefMito == "rCRS":
            # These oddball models were not found on the internet but are used in some studies deposited on ENA
            Refgenome, RefgenomeNew = ("hg37ew", "EBI37h") if  SNTypeC == "Chr" and \
                                                               SNTypeM == "chrM" else \
                                      ("hg37eht", "1K37ht") if SNTypeC == "Chr" and \
                                                               SNTypeM == "chrMT" else \
                                      ("hg37eg", "NCB37g") if  SNTypeC == "Num" and \
                                                               SNTypeM == "MT" else \
                                      (None, None)

        else:
            Refgenome = RefgenomeNew = None

        DEBUG(f'Ref Genome: {Refgenome}, by New nomenclature: {RefgenomeNew}')
        return Refgenome, RefgenomeNew

    def missing_refgenome(self, refgenome_qFN, required=True, parent=wgse.window, ref_genome=""):
        """
        Check if reference genome file exists and is correctsize; give chance to correct if not; otherwise
        report error if not corrected; return final status. Calling this routine with required=True assumes the
        ref genome is required. Otherwise simply return not missing without checking. Simplifies BAM vs CRAM checks.
        """
        missing = True      # Default return value; simply for clarity below

        if not required:        # Mainly for non-CRAMs. Just return that file is not missing if not required
            return not missing

        if self.installed_refgenome(ref_genome, refgenome_qFN):
            return not missing

        if not refgenome_qFN:           # and required (true because returned already if false)
            refgenome_qFN = self.get_refgenome_qFN(ref_genome)
            if not refgenome_qFN:       # If cannot determine file; cannot check if missing or not
                return missing

        refgenome_oFN = nativeOS(unquote(refgenome_qFN))
        refgenome_FBS = os.path.basename(refgenome_oFN)

        # Know ref genome is missing from the library now and that it is required.  So post pop-up to user and ask
        # if they want us to (1) download and install for them now, (2) give them the chance to download and install
        # while the program is paused, or (3) simply cancel the current command.  Return missing if still not there or
        # cancelled becuase if required and return missing; the caller will then cancel the command.

        # Need to reverse lookup from filename to get refgenome key; simply for error reporting
        refgenome = [key for key in self._genomes if self._genomes[key][self._final_file] == refgenome_FBS][0]
        if not refgenome:
            DEBUG(f'***Internal error: passed in {refgenome_FBS} not found in genomes.csv')

        # Give user a chance to resolve the issue by fetching the missing file (either allow us or they do it)
        message = wgse.lang.i18n["errRefGenFileMissing"]
        message = message.replace("{{FILE}}", refgenome_FBS).replace("{{REFGEN}}", refgenome)
        result = wgse_message("yesnocancel", "errRefGenFileMissingTitle", True, message)

        # User hit Cancel to decline to fix the missing refgenome. Simply return that file is missing.
        if result is None:
            return missing

        # User hit Yes to have us load the missing reference genome for them
        if result:
            get_refgenomes = nativeOS(f'{wgse.install_FP}scripts/get_refgenomes.sh')
            # Note: passing file BASE name and suffix as get_refgenome will find the URL in genomes.csv from that
            command = [wgse.bashx_oFN, get_refgenomes, refgenome_FBS, wgse.prefserver]
            run_bash_script("GenLoadGenome", command, parent=parent, direct=True)
            # The above generates a one line command to execute directly (no temp folder bash script)

        # User hit No so should have run Library.* command themselves; so hopefully should be available now ...
        # elif not result:

        # We need to check again in case the refgenome still does not exist ; also check it is processed
        if self.installed_refgenome(ref_genome, refgenome_qFN):
            return not missing

        # So an error pop-up to let user know the current command will cancel (as the file is still missing)
        message = wgse.lang.i18n["errRefGenFileStillMissing"]
        message = message.replace("{{FILE}}", refgenome_FBS).replace("{{REFGEN}}", refgenome)
        wgse_message("error", "errRefGenFileMissingTitle", True, message)

        return missing

    @staticmethod
    def refgen_liftover(ref_genome):
        """
        Return auto-selected Realign button ref. genome name and SN Count from the ref. genome name (only) passed in.
        """
        # We realign draft, early T2T style models to the T2T final
        return ("hg37",    93) if ref_genome == "hg38" else \
               ("hs37d5",  86) if ref_genome in ["hs38", "hs38a", "hs38d1", "hs38d1a", "hs38d1s",
                                                 "hs38d1v","hs38DH"] else \
               ("GRCh37", 297) if ref_genome in ["GRCh38", "GRCh38-"] else \
               ("hg38",   455) if ref_genome in ["hg37", "hg19", "hg37ew", "hg37eht", "hg37eg"] else \
               ("hs38",   195) if ref_genome in ["hs37d5", "hs37g", "hs37-", "hs37h", "hs37ht", "hs37fh", "hs37fht",
                                                 "T2Tv20", "T2Tv20a"] else \
               ("GRCh38", 639) if ref_genome in ["GRCh37", "GRCh37-"] else \
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
            return f'"{self.FP}{FBS}"'      # In reference_library; from ybrowse.org

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
        """y_only_bed_files = {
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


# The rest of this file defines a stand-alone application to process the Final Assembly Reference Genome files
#
#  A Final Assembly (FA) Human Reference Genome Processor to create keys to match to SAM/BAM/CRAM file header
#  Used to create table to figure out what FA file was used to create a BAM; that can then be used for CRAM, etc
#  FA Files are compressed in BGZF format. As we will only read, we can use standard gzip readers.
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
