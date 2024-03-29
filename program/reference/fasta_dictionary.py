from collections import OrderedDict
import typing
import pathlib
import logging

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
    def __init__(self, path: pathlib.Path) -> None:
        self._path = path
        self.header = None
        self.entries: typing.OrderedDict[str, FastaDictionaryEntry] = OrderedDict()
        self._handlers = {"@HD": self._header_process, "@SQ": self._sequence_process}
        self._load_dictionary()

    def _header_process(self, parts: typing.List[str]):
        version = None
        sorted = None
        for part in parts:
            if part.startswith("VN"):
                version = part.split(":")[1]
            elif part.startswith("SO:"):
                sorted = part.split(":")[1].strip().lower() == "sorted"
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

    def _load_dictionary(self) -> None:
        with self._path.open("rt") as f:
            for index, line in enumerate(f):
                parts = line.split()
                line_type = parts[0].strip()
                if line_type in self._handlers:
                    self._handlers[line_type](parts)
                else:
                    logging.debug(
                        f"Skipping line {index} in dictionary {self._path.name} because unrecognized type {line_type}."
                    )
