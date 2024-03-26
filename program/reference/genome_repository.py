import hashlib
import logging
from pathlib import Path
from typing import List
from .decompressor import Decompressor
from .compressor import Compressor
from .downloader import Downloader
from .file_type_checker import FileTypeChecker, Type
from .genome import Genome

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
            
    def _post_download(self, genome: Genome):
        pass
    
    def _download(self, genome: Genome):
        if genome.initial_name.exists():
            md5 = hashlib.md5(genome.initial_name.open("rb").read()).hexdigest()
            if md5 == genome.initial_md5:
                logging.info(f"{genome.code}: Already downloaded.")
                return
            logging.info(f"{genome.code}: Downloaded file is corrupted. Downloading again.")
        
        self._downloader.download(genome, None)

    def add_to_library(self, genome: Genome):
        self._get_bgzip(genome)
        self._post_download(genome)
    
    def _get_bgzip(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        if genome.final_name.exists():
            type = self._type_checker.get_type(genome.final_name)
            if type == Type.BGZIP:                
                md5 = hashlib.md5(genome.final_name.open("rb").read()).hexdigest()
                if md5 == genome.final_md5:
                    logging.info(f"{genome.code}: Already downloaded, no further conversions needed.")
                    return
                else:
                    logging.info(f"{genome.code}: Existing file is corrupted. Downloading again.")
        
        if genome.initial_name.exists():
            type = self._type_checker.get_type(genome.initial_name)
            if type == Type.BGZIP:
                md5 = hashlib.md5(genome.final_name.open("rb").read()).hexdigest()
                if md5 == genome.final_md5:
                    logging.info(f"{genome.code}: Already downloaded, no further conversions needed.")
                    return
                else:
                    logging.info(f"{genome.code}: Downloaded file is corrupted. Downloading again.")
                logging.info(f"{genome.code}: Already downloaded, no further conversions needed.")
                return
        else:
            self._download(genome)
        
        logging.info(f"{genome.code}: Converting to bgzip.")
        self.to_bgzip(genome)
