import enum
import pathlib

class Source(enum.Enum):
    AWS = 0
    NIH = 1
    NIH_ALT = 2
    EBI = 3
    EBI_ALT = 4
    UCSC = 5
    YSEQ = 6
    GOOGLE = 7
    WGSE = 8

class Genome:
    """Represent a single reference genome"""

    def __init__(self) -> None:
        self.code: str = None
        self.source: Source = None
        self.url: str = None
        self.label: str = None
        self.sn_count: str = None
        self.sn_naming: str = None
        self.description: str = None
        self.initial_name: pathlib.Path = None
        self.final_name: pathlib.Path = None