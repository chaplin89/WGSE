import enum
from .genome import Genome, Source
from pathlib import Path
import urllib.parse
from typing import Tuple, List
import csv


_STRING_TO_SOURCE = {
    "AWS": Source.AWS,
    "NIH": Source.NIH,
    "EBI": Source.EBI,
    "NIH-Alt": Source.NIH_ALT,
    "EBI-Alt": Source.EBI_ALT,
    "UCSC": Source.UCSC,
    "YSEQ": Source.YSEQ,
    "GOOG": Source.GOOGLE,
    "WGSE": Source.WGSE,
    "JHU": Source.YHU,
    "LOCAL": Source.LOCAL
}


class CsvFields(enum.IntEnum):
    INDEX = 0
    CODE = 1
    SN_COUNT = 2
    SN_NAMING = 3
    BUILD = 4
    CLASS = 5
    SOURCE = 6
    FINAL_NAME = 7
    INITIAL_NAME = 8
    URL = 9
    MENU = 10
    DESCRIPTION = 11
    FINAL_MD5 = 12
    INITIAL_MD5 = 13
    FINAL_SIZE = 14
    INITIAL_SIZE = 15


class CsvMetadataLoader:
    _EXPECTED_CSV_FIELDS = [
        "0",
        "Pyth Code (IMPORTANT: no commas in the file anywhere!)",
        "SN Count",
        "SN Name Type",
        "Build",
        "Class",
        "Source",
        "Final File Name (dot means same as initial; unique unless alt download site)",
        "Initial File Name (as downloaded; rename to final before processing)",
        "URL (unique entries; never should be a duplicate)",
        "Menu Selection",
        "Description  (If you add a column here then update zcommon.sh and referencelibrary.py)",
        "MD5Sum Final",
        "MD5Sum Initial",
        "Final Size",
        "Init Size"
    ]

    def __init__(self, csv_path: Path) -> None:
        self.metadata, self.genomes = self._load_from_csv(csv_path)

    def _load_from_csv(self, csv_path: Path) -> Tuple[List[str], List[Genome]]:
        """Load a CSV containing reference genomes information.

        Args:
            csv_path (str): Path of the CSV.

        Raises:
            RuntimeError: Raised if unable to load the CSV

        Returns:
            Tuple[list[str], list[Genome]]: A tuple containing a list of fields and a list of genome objects.
        """
        if not csv_path.exists():
            raise FileNotFoundError(f"Unable to find Genome CSV file at {csv_path}.")

        with csv_path.open() as csvfile:
            reader = csv.reader(csvfile, delimiter=",")
            rows = list(reader)

            if len(rows) < 2:
                raise RuntimeError(f"Unable to load genome reference file {csv_path}.")

            meta = rows[0]
            genomes = []

            if meta != CsvMetadataLoader._EXPECTED_CSV_FIELDS:
                raise RuntimeError(
                    f"Wrong CSV format. Expected {CsvMetadataLoader._EXPECTED_CSV_FIELDS}, got {meta}."
                )

            for row in rows[1::]:
                genome_ = Genome()
                genome_.code = row[CsvFields.CODE]
                genome_.source = _STRING_TO_SOURCE[row[CsvFields.SOURCE]]
                genome_.url = row[CsvFields.URL]
                genome_.label = row[CsvFields.MENU]
                genome_.sn_count = row[CsvFields.SN_COUNT]
                genome_.sn_naming = row[CsvFields.SN_NAMING]
                genome_.description = row[CsvFields.DESCRIPTION]
                genome_.final_name = ""  # row[CsvFields.FINAL_NAME]
                genome_.initial_name = ""  # Path(row[CsvFields.INITIAL_NAME])
                genome_.initial_size = int(row[CsvFields.INITIAL_SIZE])
                genome_.final_size = int(row[CsvFields.FINAL_SIZE])
                genome_.initial_md5 = row[CsvFields.INITIAL_MD5]
                genome_.final_md5 = row[CsvFields.FINAL_MD5]
                
                if genome_.initial_name == "":
                    extension = self.get_compression_extension(genome_.url)
                    genome_.initial_name = Path(
                        f"{genome_.code}_{genome_.source.name}.fa{extension}"
                    )
                if genome_.final_name == "":
                    genome_.final_name = self.determine_target_name(
                        genome_.initial_name
                    )

                genomes.append(genome_)
            return meta, genomes

    def determine_target_name(self, file: Path):
        if file.suffix in [".bz2", ".bz", ".zip", ".7z"]:
            gz_compressed = file.with_suffix(".gz")
        elif file.suffix != ".gz":
            # Likely uncompressed, just add .gz
            gz_compressed = Path(str(file) + ".gz")
        else:
            # gzip case, same name
            gz_compressed = file
        
        return gz_compressed

    def get_compression_extension(self, url: str):
        parsed = urllib.parse.urlparse(url)
        recognized_extensions = [
            ".zip",
            ".7z",
            ".gz",
            ".bz2"
        ]
        
        for extension in recognized_extensions:
            if parsed.path.endswith(extension):
                return extension
        return ""

    def filter(
        self, id: str = None, source: Source = None, file_name: str = None
    ) -> List[Genome]:
        """Return a list of genome objects filtered by parameters.

        Args:
            id (str, optional): ID of the genome. Defaults to None.
            source (Source, optional): Source of the genome. Defaults to None.
            file_name (str, optional): File name of the genome. Defaults to None.

        Returns:
            List[Genome]: List of genome objects matching the criteria.
        """
        filtered = self.genomes.copy()
        if id is not None:
            filtered = [x for x in filtered if x.code == id]
        if source is not None:
            filtered = [x for x in filtered if x.source == source]
        if file_name is not None:
            filtered = [x for x in filtered if x.initial_name.name == file_name]
        return filtered

    def single(
        self, id: str = None, source: Source = None, file_name: str = None
    ) -> Genome:
        """Same as filter but raises an exception if not exactly a single genome is returned.

        Args:
            id (str, optional): ID of the genome. Defaults to None.
            source (Source, optional): Source of the genome. Defaults to None.
            file_name (str, optional): File name of the genome. Defaults to None.

        Raises:
            RuntimeError: If after the filtering not exactly a single genome is returned.

        Returns:
            Genome: The specific genome requested by the caller.
        """
        filtered = self.filter(id, source, file_name)
        if len(filtered) != 1:
            raise RuntimeError(f"Expected one genome, found {len(filtered)}")
        return filtered
