from .external import External
from pathlib import Path


class Compressor:
    def __init__(self, external: External) -> None:
        self._external = External

    def compress(self, input_file: Path) -> Path:
            return self._external.bgzip(str(input_file))