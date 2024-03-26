from .samtools import Samtools
from pathlib import Path


class Compressor:
    def __init__(self, samtools: Samtools) -> None:
        self._samtools = samtools

    def compress(self, input_file: Path) -> Path:
            return self._samtools.bgzip_compress(str(input_file))