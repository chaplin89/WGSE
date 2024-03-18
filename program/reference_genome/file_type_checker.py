from pathlib import Path
from samtools import Samtools
import logging
import enum


class Type(enum.Enum):
    ZIP = 0
    BZIP = 1
    SEVENZIP = 2
    GZIP = 3
    BGZIP = 4
    DECOMPRESSED = 5


class FileTypeChecker:
    _EXT_TO_TYPE = {
        ".7z": Type.SEVENZIP,
        ".zip": Type.ZIP,
        ".bz": Type.BZIP,
        ".bz2": Type.BZIP,
        ".gz": Type.GZIP,
        ".fasta": Type.DECOMPRESSED,
        ".fna": Type.DECOMPRESSED,
        ".fa": Type.DECOMPRESSED,
    }

    _HTSFILE_TO_TYPE = {
        "BGZF-compressed": Type.BGZIP,
        "gzip-compressed": Type.GZIP,
        "RAZF-compressed": Type.GZIP,
        "FASTA sequence": Type.DECOMPRESSED,
    }

    def __init__(self, samtools: Samtools) -> None:
        self._samtools = samtools

    def get_type(self, file: Path) -> Type | None:
        # Extensions can be wrong or misleading; Use htsfile and
        # eventually fallback on extension.

        file_type = self._samtools.get_file_type(file)
        logging.debug(f"htsfile returned {file_type}.")
        for key in FileTypeChecker._HTSFILE_TO_TYPE.keys():
            if key in file_type:
                return FileTypeChecker._HTSFILE_TO_TYPE[key]

        # TODO: starting from .gz in _EXT_TO_TYPE, htsfile should really
        # be able to give something useful. If it isn't the case, it could
        # mean the file is corrupted. Unsure what's the best logic here.
        extension = file.suffix
        if extension in FileTypeChecker._EXT_TO_TYPE:
            return FileTypeChecker._EXT_TO_TYPE[extension]
        return None

    def need_decompression(self, file: Path):
        type = self.get_type(file)
        return type != Type.BGZIP and type != Type.DECOMPRESSED

    def need_compression(self, file: Path):
        type = self.get_type(file)
        return type != Type.BGZIP
