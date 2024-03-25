import logging
from pathlib import Path
from typing import List
from decompressor import Decompressor
from compressor import Compressor
from downloader import Downloader
from file_type_checker import FileTypeChecker, Type
from genome import Genome

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
    ) -> None:
        self._reference_genomes = reference_genomes
        self._library_path = library_path
        self._type_checker = type_checker
        self._downloader = downloader_
        self._decompressor = decompressor_
        self._compressor = compressor_
        self._rebase_reference_genomes()

    def _rebase_reference_genomes(self):
        for genome in self._reference_genomes:
            genome.initial_name = self._library_path.joinpath(genome.initial_name)
            genome.final_name = self._library_path.joinpath(genome.final_name)

    def to_bgzip(self, genome: Genome):
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
                logging.info(f"{genome.code}: bgzip file found. No further conversions needed.")
                return
        
        if type != Type.DECOMPRESSED:
            logging.info(f"{genome.code}: Decompressing.")
            decompressed = self._decompressor.decompress(genome.initial_name)
        else:
            logging.info(f"{genome.code}: Is already decompressed. Skipping decompression.")
            decompressed = genome.initial_name

        logging.info(f"{genome.code}: Compressing to bgzip.")
        compressed = self._compressor.compress(decompressed)
        # decompressed no longer exist at this point
        
        if genome.initial_name.exists() and compressed != genome.initial_name:
            genome.initial_name.unlink()

        if compressed != genome.final_name:
            compressed.rename(genome.final_name)

    def add_to_library(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        logging.info(f"{genome.code}: Adding to library.")
        if self._downloader.need_download(genome):
            logging.info(f"{genome.code}: Downloading Genome.")
            self._downloader.download(genome, None)
        else:
            logging.info(f"{genome.code}: Already downloaded.")

        type = self._type_checker.get_type(genome.initial_name)
        logging.info(f"{genome.code}: Determined type {type.name}.")

        if type == Type.BGZIP:
            logging.info(f"{genome.code}: Already in bgzip format, no further conversion required.")
        else:
            try:
                logging.info(f"{genome.code}: Converting to bgzip.")
                self.to_bgzip(genome)
            except Exception as e:
                logging.error(e)
                raise

        # t = samtools.fasta_index(genome.final_name)
        # t = samtools.make_dictionary(genome.final_name)
