import enum
import hashlib
import logging
from pathlib import Path
from typing import List

from .fasta_file import FastaFile
from .external import GzipAction, External
from .decompressor import Decompressor
from .compressor import Compressor
from .downloader import Downloader
from .file_type_checker import FileTypeChecker, Type
from .genome import Genome
from .n_stats_files import NStatsFiles
from .csv_metadata_loader import CsvMetadataLoader

class FileStatus(enum.Enum):
    Available = 0,
    NotAvailable = 1,
    Corrupt = 2,
    Valid = 3,

class GenomeRepository:
    """Manage a repository of reference genome files"""

    def build(csv: Path, root: Path, external: Path = None):
        external = External(external)
        type_checker = FileTypeChecker(external)
        return GenomeRepository(
            CsvMetadataLoader(csv).genomes,
            root,
            type_checker,
            Downloader(),
            Compressor(external),
            Decompressor(type_checker, external),
            external
        )

    def __init__(
        self,
        reference_genomes: List[Genome],
        library_path: Path,
        type_checker: FileTypeChecker,
        downloader_: Downloader,
        compressor_: Compressor,
        decompressor_: Decompressor,
        external_: External,
    ) -> None:
        self._reference_genomes = reference_genomes
        self._library_path = library_path
        self._type_checker = type_checker
        self._downloader = downloader_
        self._decompressor = decompressor_
        self._compressor = compressor_
        self._external = external_
        self._rebase_reference_genomes()

    def _rebase_reference_genomes(self):
        for genome in self._reference_genomes:
            genome.initial_name = self._library_path.joinpath(genome.initial_name)
            genome.final_name = self._library_path.joinpath(genome.final_name)

    def _to_bgzip(self, genome: Genome):
        type = self._type_checker.get_type(genome.initial_name)

        if type == Type.BGZIP:
            return


        if type != Type.DECOMPRESSED:
            logging.info(f"{genome.code}-{genome.source.name}: Decompressing.")
            self._decompressor.decompress(genome.initial_name, genome.decompressed)
        else:
            logging.info(
                f"{genome.code}-{genome.source.name}: Is already decompressed. Skipping decompression."
            )

        logging.info(f"{genome.code}-{genome.source.name}: Compressing to bgzip.")
        self._compressor.compress(genome.decompressed, genome.final_name)
        # Decompressed no longer exist at this point

        # Need to consider that initial and final name may be identical.
        if genome.initial_name.exists() and genome.final_name != genome.initial_name:
            genome.initial_name.unlink()

    def filter(self, id=None, source=None):
        filtered = self._reference_genomes
        if id is not None:
            filtered = [x for x in filtered if x.code.startswith(id)]
        if source is not None:
            filtered = [x for x in filtered if x.source.name == source]
        return filtered

    def _post_download(self, genome: Genome):
        logging.info(f"{genome.code}-{genome.source.name}: Starting post-download tasks.")
        if not genome.gzi.exists():
            logging.info(f"{genome.code}-{genome.source.name}: Generating bgzip index.")
            index = self._external.bgzip(genome.final_name, GzipAction.Reindex)
            if index != genome.gzi:
                index.rename(genome.gzi)
        else:
            logging.info(f"{genome.code}-{genome.source.name}: bgzip index exists.")
        if not genome.dict.exists():
            logging.info(f"{genome.code}-{genome.source.name}: Generating .dict file.")
            self._external.make_dictionary(genome.final_name, genome.dict)
        else:
            logging.info(f"{genome.code}-{genome.source.name}: .dict file exists.")
        if not all([genome.bed.exists(), genome.nbin.exists(), genome.nbuc.exists()]):
            logging.info(f"{genome.code}-{genome.source.name}: Generating Ns stats files.")
            fasta_file = FastaFile(genome)
            ub = NStatsFiles(fasta_file)
            ub.generate_stats()
        else:
            logging.info(f"{genome.code}-{genome.source.name}: Ns stats files exist.")
        logging.info(f"{genome.code}-{genome.source.name}: Post-download tasks ended.")

    def _check_file_md5(self, path: Path, md5: str):
        if not path.exists():
            raise FileNotFoundError(f"Unable to find file {str(path)} to calculate MD5.")
        md5_algorithm = hashlib.md5()
        with path.open("rb") as f:
            while chunk := f.read(4096):
                md5_algorithm.update(chunk)

        if md5_algorithm.hexdigest() == md5:
            return True
        return False

    def add(self, genome: Genome):
        self._get_bgzip(genome)
        self._post_download(genome)
    
    def delete(self, genome:Genome) -> List[Path]:
        deleted = []
        for file in genome.all:
            if file.exists():
                file.unlink()
                deleted.append(file)
        return deleted
    
    def _check_file(self, file: Path, md5: str):
        if file.exists():
            self._check_file_md5(file, md5)
        if file.exists():
            md5_matches = self._check_file_md5(file, md5)
            if md5_matches:
                return FileStatus.Valid
            else:
                return FileStatus.Corrupt
        return FileStatus.NotAvailable
    
    def check(self, genome:Genome):
        status = {}
        final = self._check_file(genome.final_name, genome.final_md5)
        initial = self._check_file(genome.initial_name, genome.initial_md5)
        if genome.initial_name == genome.final_name:
            if final == FileStatus.Valid or initial == FileStatus.Valid:
                status[genome.initial_name] = FileStatus.Valid
            else:
                status[genome.initial_name] = final
        else:
            status[genome.initial_name] = initial
            status[genome.final_name] = final
        
        for file in genome.all:
            if file not in status:
                if file.exists():
                    status[file] = FileStatus.Available
                else:
                    status[file] = FileStatus.NotAvailable
        return status
        

    def _get_bgzip(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        if genome.final_name.exists():
            type = self._type_checker.get_type(genome.final_name)
            if type == Type.BGZIP:
                logging.info(f"{genome.code}-{genome.source.name}: found a file in target file format. Checking if the file is corrupt.")
                md5_matches = self._check_file_md5(genome.final_name, genome.final_md5)
                if md5_matches:
                    logging.info(f"{genome.code}-{genome.source.name}: File is OK. Nothing left to do.")
                    return
                else:
                    logging.info(f"{genome.code}-{genome.source.name}: File is corrupt. Trying to add again to the repository.")

        need_download = True
        if genome.initial_name.exists():
            logging.info(f"{genome.code}-{genome.source.name}: found a file in initial file format. Checking MD5.")
            md5_matches = self._check_file_md5(genome.initial_name, genome.initial_md5)
            if md5_matches:
                logging.info(f"{genome.code}-{genome.source.name}: File is OK. No need to download it again.")
                need_download = False
            elif md5_matches == False:
                logging.info(
                    f"{genome.code}-{genome.source.name}: {genome.initial_name.name} is corrupted. Downloading again."
                )
        else:
            logging.info(f"{genome.code}-{genome.source.name}: File not found. Downloading.")

        if need_download:
            self._downloader.download(genome, False)
        self._to_bgzip(genome)