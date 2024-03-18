import enum
from pathlib import Path
from typing import List

class Source(enum.Enum):
    AWS = 0
    NIH = 1
    EBI = 2
    NIH_ALT = 3
    EBI_ALT = 4
    UCSC = 5
    YSEQ = 6
    GOOG = 7
    WGSE = 8

class Genome:
    """Represent a single reference genome"""

    def __init__(self) -> None:
        self.code: str = None
        self.source: Source = None
        self.final_name: Path = None
        self.initial_name: Path = None
        self.uncompressed: Path = None
        self.bgzip_compressed: Path = None
        self.url: str = None
        self.label: str = None
        self.sn_count: str = None
        self.sn_naming: str = None
        self.description: str = None