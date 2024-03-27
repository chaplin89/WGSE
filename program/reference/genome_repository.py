import hashlib
import logging
from pathlib import Path
from typing import List

from .fasta_file import FastaFile
from .samtools import BgzipAction, Samtools
from .decompressor import Decompressor
from .compressor import Compressor
from .downloader import Downloader
from .file_type_checker import FileTypeChecker, Type
from .genome import Genome
from .n_statistics_files import NStatisticsFiles


class GenomeRepository:
    """Manage a repository of reference genome files"""

    def __init__(
        self,
        reference_genomes: List[Genome],
        library_path: Path,
        type_checker: FileTypeChecker,
        downloader_: Downloader,
        compressor_: Compressor,
        decompressor_: Decompressor,
        samtools_ : Samtools,
    ) -> None:
        self._reference_genomes = reference_genomes
        self._library_path = library_path
        self._type_checker = type_checker
        self._downloader = downloader_
        self._decompressor = decompressor_
        self._compressor = compressor_
        self._samtools = samtools_
        self._rebase_reference_genomes()

    def _rebase_reference_genomes(self):
        for genome in self._reference_genomes:
            genome.initial_name = self._library_path.joinpath(genome.initial_name)
            genome.final_name = self._library_path.joinpath(genome.final_name)

    def _to_bgzip(self, genome: Genome):
        type = self._type_checker.get_type(genome.initial_name)

        if type == Type.RAZF_GZIP:
            logging.info(f"{genome.code}: RAZF is not supported. Skipping.")
            # TODO: figure out how to decompress this format.
            return

        if type == Type.BGZIP:
            return

        if genome.final_name.exists():
            final_type = self._type_checker.get_type(genome.final_name)
            if final_type == Type.BGZIP:
                logging.info(
                    f"{genome.code}: bgzip file found. No further conversions needed."
                )
                return

        if type != Type.DECOMPRESSED:
            logging.info(f"{genome.code}: Decompressing.")
            decompressed = self._decompressor.decompress(genome.initial_name)
        else:
            logging.info(
                f"{genome.code}: Is already decompressed. Skipping decompression."
            )
            decompressed = genome.initial_name

        logging.info(f"{genome.code}: Compressing to bgzip.")
        compressed = self._compressor.compress(decompressed)
        # Decompressed no longer exist at this point

        # Need to consider that initial and final name may be identical.
        if genome.initial_name.exists() and compressed != genome.initial_name:
            genome.initial_name.unlink()

        if compressed != genome.final_name:
            compressed.rename(genome.final_name)

    def _post_download(self, genome: Genome):
        if not genome.gzi.exists():
            index = self._samtools.bgzip(genome.final_name, BgzipAction.Reindex)
            if index != genome.gzi:
                index.rename(genome.gzi)
        if not genome.dict.exists():
            self._samtools.make_dictionary(genome.final_name, genome.dict)
        if not all([genome.bed.exists(), genome.nbin.exists(), genome.nbuc.exists()]):
            fasta_file = FastaFile(genome)
            ub = NStatisticsFiles(fasta_file)
            ub.generate_stats()

    def _check_file_md5(self, path: Path, md5: str):
        if path.exists():
            md5_algorithm = hashlib.md5()
            with path.open("rb") as f:
                while chunk := f.read(4096):
                    md5_algorithm.update(chunk)

            if md5_algorithm.hexdigest() == md5:
                return True
            return False
        return None

    def _download(self, genome: Genome):
        self._downloader.download(genome, None)

    def add_to_library(self, genome: Genome):
        self._get_bgzip(genome)
        self._post_download(genome)

    def _get_bgzip(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        md5_matches = self._check_file_md5(genome.final_name, genome.final_md5)

        if md5_matches:
            type = self._type_checker.get_type(genome.final_name)
            if type == Type.BGZIP:
                logging.info(f"{genome.code}: file downloaded, no conversion needed.")
                return

        need_download=True
        md5_matches = self._check_file_md5(genome.initial_name, genome.initial_md5)
        if md5_matches:
            type = self._type_checker.get_type(genome.initial_name)
            if type == Type.BGZIP:
                # This should never happen as we have a bgzip file we should be falling in the previous if
                raise RuntimeError(
                    f"{genome.code}: {genome.final_name.name} not present but {genome.initial_name.name} in bgzip format."
                )
            logging.info(f"{genome.code}: Already downloaded.")
            need_download = False
        elif md5_matches == False:
            logging.info(
                f"{genome.code}: {genome.initial_name.name} is corrupted. Downloading again."
            )
        elif md5_matches == None:
            logging.info(
                f"{genome.code}: File not present. Downloading."
            )
            
        if need_download:
            self._download(genome)        
        self._to_bgzip(genome)
