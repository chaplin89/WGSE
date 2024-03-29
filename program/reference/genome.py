import enum
import pathlib
import typing


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
    LOCAL = 10


class Genome:
    """Represent a single reference genome"""

    def __init__(self, path: pathlib.Path = None) -> None:
        self.code: str = None
        self.source: Source = None
        self.url: str = None
        self.label: str = None
        self.sn_count: str = None
        self.sn_naming: str = None
        self.description: str = None
        self.initial_name: pathlib.Path = None
        self.final_name: pathlib.Path = path
        self.initial_md5: str = None
        self.final_md5: str = None
        self.initial_size: int = None
        self.final_size: int = None

    @property
    def no_exts(self) -> str:
        if self.final_name is None:
            return None
        name = str(self.final_name)
        extensions = "".join(self.final_name.suffixes)
        return name.rstrip(extensions)
    
    @property
    def all(self) -> typing.List[pathlib.Path]:
        return [
            self.initial_name,
            self.final_name,
            self.decompressed,
            self.gzi,
            self.nbin,
            self.nbuc,
            self.bed,
            self.fai,
            self.dict
        ]

    @property
    def decompressed(self) -> pathlib.Path:
        if self.final_name is None:
            return None
        return pathlib.Path(self.no_exts + ".fa")
    
    @property
    def gzi(self):
        if self.final_name is None:
            return None
        return pathlib.Path(str(self.final_name) + ".gzi")

    @property
    def fai(self) -> pathlib.Path:
        if self.final_name is None:
            return None
        return pathlib.Path(str(self.final_name) + ".fai")
    
    @property
    def nbin(self) -> pathlib.Path:
        if self.final_name is None:
            return None
        return pathlib.Path(self.no_exts + "_nbin.csv")

    @property
    def nbuc(self) -> pathlib.Path:
        if self.final_name is None:
            return None
        return pathlib.Path(self.no_exts + "_nbuc.csv")

    @property
    def bed(self) -> pathlib.Path:
        if self.final_name is not None:
            return pathlib.Path(self.no_exts + "_nreg.bed")
        return None

    @property
    def dict(self) -> pathlib.Path:
        if self.final_name is None:
            return None
        return pathlib.Path(self.no_exts + ".dict")

    def __repr__(self) -> str:
        return f"{self.code} from {self.source.name}"
    
    def __str__(self) -> str:
        return self.__repr__()