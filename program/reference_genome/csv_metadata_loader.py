
import enum
from genome import Genome, Source
from pathlib import Path
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
}

class CsvFields(enum.IntEnum):
    CODE = 0
    SOURCE = 1
    FINAL_NAME = 2
    INITIAL_NAME = 3
    URL = 4
    LABEL = 5
    SN_COUNT = 6
    SN_NAMING = 7
    DESCRIPTION = 8

class CsvMetadataLoader:
    _EXPECTED_CSV_FIELDS = [
        "Pyth Code",
        "Source",
        "Final File Name",
        "Initial File Name",
        "URL",
        "Library Menu Label",
        "SN Count",
        "SN Naming",
        "Description",
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
                genome_.label = row[CsvFields.LABEL]
                genome_.sn_count = row[CsvFields.SN_COUNT]
                genome_.sn_naming = row[CsvFields.SN_NAMING]
                genome_.description = row[CsvFields.DESCRIPTION]
                genome_.final_name = row[CsvFields.FINAL_NAME]
                genome_.initial_name = Path(row[CsvFields.INITIAL_NAME])
                genome_.file = Path(genome_.code + "_" + genome_.source.name + "".join(genome_.initial_name.suffixes))
                
                genomes.append(genome_)
            return meta, genomes

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