from .external import External
from pathlib import Path


class Compressor:
    def __init__(self, external: External) -> None:
        self._external = external

    def compress(self, input_file: Path, output_file: Path) -> Path:
            return self._external.bgzip(input_file, output_file)