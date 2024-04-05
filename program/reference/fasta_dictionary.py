from collections import OrderedDict
import enum
import typing
import pathlib
import logging


class MTDNANameType(enum.Enum):
    chrM = 0
    MT = 1
    chrMT = 2
    M = 3
    Accession = 4


class ChromosomeNameType(enum.Enum):
    Chr = 0
    Number = 1
    Accession = 2


class Gender(enum.Enum):
    Male = 0
    Female = 1


class Sorting(enum.Enum):
    Coordinate = 0
    Name = 1
    Unsorted = 2


class FastaDictionaryEntry:
    def __init__(self, name: str, length: int, md5: str, uri: str) -> None:
        self.name: str = name
        self.length: int = length
        self.md5: str = md5
        self.uri: str = uri


class FastaDictionaryHeader:
    def __init__(self, version: str, sorted_: bool) -> None:
        self.version: str = version
        self.sorted: bool = sorted_


class FastaDictionary:
    def __init__(self, lines) -> None:
        self.header = None
        self.entries: typing.OrderedDict[str, FastaDictionaryEntry] = OrderedDict()
        self._handlers = {"@HD": self._header_process, "@SQ": self._sequence_process}
        self._load(lines)


    def _header_process(self, parts: typing.List[str]):
        version = None
        sorted = None
        for part in parts:
            if part.startswith("VN"):
                version = part.split(":")[1]
            elif part.startswith("SO:"):
                sorted = self._is_sorted(part)
        self.header = FastaDictionaryHeader(version, sorted)

    def _sequence_process(self, parts: typing.List[str]):
        name = None
        length = None
        md5 = None
        uri = None

        for part in parts:
            if part.startswith("SN"):
                name = part[3::].strip()
            elif part.startswith("LN"):
                length = int(part[3::].strip())
            elif part.startswith("M5"):
                md5 = part[3::].strip()
            elif part.startswith("UR"):
                uri = part[3::].strip()
        if name == None:
            raise RuntimeError("Unable to find the name of the sequence.")
        if length == None:
            raise RuntimeError("Unable to find the length of the sequence.")
        entry = FastaDictionaryEntry(name, length, md5, uri)
        self.entries[name] = entry

    def _load(self, lines):
        for index, line in enumerate(lines):
            line = line.strip()
            if line == "":
                continue
            parts = line.split()
            line_type = parts[0].strip()
            if line_type in self._handlers:
                self._handlers[line_type](parts)
            else:
                logging.debug(
                    f"Skipping line {index} in dictionary because unrecognized type {line_type}."
                )

    def load_from_file(path: pathlib.Path) -> 'FastaDictionary':
        if path is None:
            raise ValueError("path cannot be None.")
        if not path.exists():
            raise FileNotFoundError(f"Unable to find file: {str(path)}")

        with path.open("rt") as f:
            return FastaDictionary(f)

    def _is_sorted(self, value: str) -> Sorting:
        value = value.split(":")[1].strip().lower()
        sorting = None
        if "coordinate" in value:
            sorting = Sorting.Coordinate
        elif "unsorted" in value:
            sorting = Sorting.Name
        elif "unknown" in value:
            sorting = Sorting.Unsorted
        else:
            raise RuntimeError(f"Unable to get sorting for {sorting}")
        return sorting

    def chromosome_name_type(self):
        patterns_type = {
            ChromosomeNameType.Accession: ["CM0", "CP", "J0", "NC_"],
            ChromosomeNameType.Chr: ["chr"],
            ChromosomeNameType.Number: ["1"],
        }

        for type, patterns in patterns_type.items():
            if any([x in self.entries.keys() for x in patterns]):
                return type
        return None

    def mtdna_name_type(self):
        patterns_type = {
            MTDNANameType.chrMT: ["chrMT"],
            MTDNANameType.MT: ["MT"],
            MTDNANameType.chrM: ["chrM"],
            MTDNANameType.M: ["M"],
        }

        for type, patterns in patterns_type.items():
            if any([x in self.entries.keys() for x in patterns]):
                return type

        if self.chromosome_name_type() == ChromosomeNameType.Accession:
            return MTDNANameType.Accession
        return None

    def sequence_count(self):
        return len(self.entries)
