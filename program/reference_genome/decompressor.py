import logging
import zipfile
import gzip
import bz2
import shutil
from genome import Genome
from samtools import Samtools
from pathlib import Path


class Decompressor:
    def __init__(self, samtools: Samtools) -> None:
        self._samtools = samtools

        self.handlers = {
            "gz": Decompressor.gz,
            "zip": Decompressor.zip,
            "7z": Decompressor.sevenzip,
            "bz": Decompressor.bzip,
            "bz2": Decompressor.bzip,
        }

        self.ignored_format = ["fa", "fna", "fasta"]

    def gz(self, file: Path):
        logging.debug(f"Decompressing file {file}. gzip compression detected.")
        target_file = file.name.strip(file.suffix)
        target_file = file.parent.joinpath(target_file)

        with gzip.open(file, "rb") as f_in:
            with open(target_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def sevenzip(self, file: Path):
        logging.debug(f"Decompressing file {file}. 7z compression detected.")
        raise NotImplementedError()

    def bzip(self, file: Path):
        logging.debug(f"Decompressing file {file}. bzip compression detected.")
        raise NotImplementedError()

    def zip(self, file: Path):
        logging.debug(f"Decompressing file {file}. zip compression detected.")
        with zipfile.ZipFile(file, "r") as f:
            f.extractall(file.parent)

    def need_decompression(self, genome: Genome) -> bool:
        """Check if the genome already has a file with target
          compression format (bgzip) or if it is uncompressed.

        Args:
            genome (Genome): Genome file to check.

        Raises:
            FileNotFoundError: File to check does not exist.

        Returns:
            bool: True if we don't have a file with target 
            compression or uncompressed. False otherwise.
        """
        if not genome.initial_name.exists():
            raise FileNotFoundError(
                f"Unable to find reference genome file {genome.initial_name}."
            )

        if genome.uncompressed.exists():
            return False
        
        if self._is_bgzip_compressed(genome.initial_name):
            return False
        
        if self._is_uncompressed(genome.initial_name):
            return False
        return True

    def _is_bgzip_compressed(self, file: Path):
        """Uses samtools to check if a file is bgzip compressed.
        
        Args:
            file (Path): File to check

        Returns:
            bool: True if bgzip compressed, False otherwise.
        """
        # Usually bgzip files have .gz extension but they can be 
        # also .fasta, .fna, .fna. This is why this check is done 
        # with samtools and not by looking at the extension.
        file_type = self._samtools.get_file_type(file)

        if "FASTA BGZF" in file_type:
            return True
        return False

    def _is_uncompressed(self, file: Path):
        """Check if the extension of a file indicates an uncompressed format.

        Args:
            file (Path): file to check

        Returns:
            bool: True if uncompressed, false otherwise.
        """
        # This method rely on checking the extension and it makes sense only
        # if the file was already determined to not be bgzip compressed.
        if file.suffix in self.ignored_format:
            return True
        return False

    def decompress(self, genome: Genome):
        """Decompress a genome file to genome.uncompressed file path.

        Args:
            genome (Genome): Genome to decompress.

        Raises:
            FileNotFoundError: Input compressed genome file is not found.
            RuntimeError: Extension is not recognized.
        """
        if not genome.initial_name.exists():
            raise FileNotFoundError(
                f"Unable to find reference genome file {genome.initial_name} to decompress."
            )
        extension = genome.initial_name.suffix
        extension = extension.lstrip(".")
        if extension not in self.handlers:
            raise RuntimeError(
                f"Trying to decompress a reference genome file with an unknown file extension: {genome.initial_name}"
            )
        
        handler = self.handlers[extension]
        handler(self, genome.initial_name)
