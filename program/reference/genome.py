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
    YHU = 9

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
        self.gzi: pathlib.Path = None
        self.nbuc : pathlib.Path = None
        self.nreg : pathlib.Path = None
        self.bed : pathlib.Path = None
        self.nbin : pathlib.Path = None
        self.fai : pathlib.Path = None
        self.dict : pathlib.Path = None
        self.initial_md5: str = None
        self.final_md5: str = None
        self.initial_size: int = None
        self.final_size: int = None
