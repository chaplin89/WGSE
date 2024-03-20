import logging
import zipfile
import gzip
import bz2
import shutil
from genome import Genome
from samtools import Samtools
from pathlib import Path
from file_type_checker import Type, FileTypeChecker


class Decompressor:
    def __init__(self, type_checker: FileTypeChecker) -> None:
        self._type_checker = type_checker

        self.handlers = {
            Type.GZIP: Decompressor.gz,
            Type.ZIP: Decompressor.zip,
            Type.SEVENZIP: Decompressor.sevenzip,
            Type.BZIP: Decompressor.bzip,
        }

    def gz(self, file: Path):
        logging.debug(f"Decompressing file {file}. gzip compression detected.")
        target_file = self.determine_target_name(file)

        with gzip.open(str(file), "rb") as f_in:
            with open(target_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return target_file

    def sevenzip(self, file: Path):
        logging.debug(f"Decompressing file {file}. 7z compression detected.")
        raise NotImplementedError()

    def bzip(self, file: Path):
        logging.debug(f"Decompressing file {file}. bzip compression detected.")
        target_file = self.determine_target_name(file)
        raise NotImplementedError()
        return target_file

    def zip(self, file: Path):
        logging.debug(f"Decompressing file {file}. zip compression detected.")
        target_file = self.determine_target_name(file)
        
        raise NotImplementedError()
        with zipfile.ZipFile(file, "r") as f:
            f.extractall(file.parent)
        return target_file

    def determine_target_name(self, file: Path):
        if file.suffix in [".bz2", ".bz", ".zip", ".7z", ".gz"]:
            return Path(file.parent, file.stem)
        raise RuntimeError(f"Unable to determine an appropriate decompressed filename for {file}")

    def decompress(self, genome: Genome):
        """Decompress a genome file.

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
        type = self._type_checker.get_type(genome.initial_name)
        if type not in self.handlers:
            raise RuntimeError(
                f"Trying to decompress a reference genome file with an unknown file extension: {genome.initial_name}"
            )
        
        handler = self.handlers[type]
        return handler(self, genome.initial_name)
