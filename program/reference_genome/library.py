import csv
import enum
import os
from pathlib import Path
from typing import List, Union, Tuple
from decompressor import Decompressor
from compressor import Compressor
from downloader import Downloader
from genome import Genome, Source
from samtools import Samtools
import argparse

_STRING_TO_SOURCE = {
    "AWS": Source.AWS,
    "NIH": Source.NIH,
    "EBI": Source.EBI,
    "NIH-Alt": Source.NIH_ALT,
    "EBI-Alt": Source.EBI_ALT,
    "UCSC": Source.UCSC,
    "YSEQ": Source.YSEQ,
    "GOOG": Source.GOOG,
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


class Library:
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

    def __init__(self, csv_path: str, library_path: str) -> None:
        """Initialize the object and load the CSV.

        Args:
            csv_path (str): Path of the CSV
            library_path (str): Path to the library

        Raises:
            RuntimeError: If the CSV format is not matching a known format.
        """
        self.library_path = Path(library_path)
        self._meta, self._genome = self._load_csv(csv_path)

    def _load_csv(self, csv_path) -> Tuple[list[str], list[Genome]]:
        """Load a CSV containing reference genomes information.

        Args:
            genome_csv (str): Path of the CSV.

        Raises:
            RuntimeError: Raised if unable to load the CSV

        Returns:
            Tuple[list[str], list[Genome]]: A tuple containing reference genome info.
        """
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Unable to find Genome CSV file at {csv_path}.")

        with open(csv_path) as csvfile:
            reader = csv.reader(csvfile, delimiter=",")
            rows = list(reader)

            if len(rows) < 2:
                raise RuntimeError(f"Unable to load genome reference file {csv_path}.")

            meta = rows[0]
            genomes = []

            if meta != Library._EXPECTED_CSV_FIELDS:
                raise RuntimeError(
                    f"Wrong CSV format. Expected {Library._EXPECTED_CSV_FIELDS}, got {self._meta}."
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
                genome_.initial_name = self.library_path.joinpath(row[CsvFields.INITIAL_NAME])
                genome_.final_name = self.library_path.joinpath(row[CsvFields.FINAL_NAME])
                genome_.uncompressed = genome_.initial_name
                if genome_.initial_name.suffix in [".zip", ".7z", ".gz", ".bz", "bz2"]:
                    genome_.uncompressed=Path(str(genome_.initial_name).rstrip(genome_.initial_name.suffix))
                genome_.bgzip_compressed = None
                genomes.append(genome_)
            return meta, genomes

    def filter(
        self, id=None, source=None, initial_name=None, final_name=None
    ) -> Union[List[Genome] | Genome]:
        """Return a list of Genome objects filtered by arguments.

        Args:
            id (str, optional): ID of the genome. Defaults to None.
            source (str, optional): Source of the genome. Defaults to None.
            initial_name (str, optional): Source file name. Defaults to None.
            final_name (str, optional): Destination file name. Defaults to None.

        Returns:
            Union[List[Genome] | Genome]: A filtered list or a single Genome objects.
        """
        filtered = self._genome.copy()
        if id is not None:
            filtered = [x for x in filtered if x.code == id]
        if source is not None:
            filtered = [x for x in filtered if x.source == source]
        if initial_name is not None:
            filtered = [x for x in filtered if x.initial_name == initial_name]
        if final_name is not None:
            filtered = [x for x in filtered if x.final_name == final_name]

        if len(filtered) == 1:
            return filtered[0]
        return filtered

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manage a library of reference genomes.', 
                                     prog=Path(__file__).name,
                                     usage='%(prog)s [options]')
    
    parser.add_argument('--csv', type=Path, metavar="path",
                        help='A CSV file containing a description on known reference genomes.',
                        default=Path("reference", "seed_genomes.csv"))
    
    parser.add_argument('--root', type=Path, metavar="path",
                        help='Path to the root folder of the library.',
                        default=Path("reference", "genomes"))
    
    parser.add_argument('--samtool', type=str, metavar="path",
                        help='Path to the folder containing samtool executables.',
                        default=Path("cygwin64","usr","local","bin"))

    download_subparser = parser.add_subparsers(help='Download a reference genome file')
    download_action = download_subparser.add_parser('download', help='Download a reference genome file')
    download_action.add_argument('id', type=str, help='ID of the reference genome to download.')

    #list_subparser = parser.add_subparsers(help='list help')
    #list_action = list_subparser.add_parser('list', help='List known reference genomes.')

    args = parser.parse_args()
    parser.print_help()

    library = Library(args.csv, args.root)
    samtools = Samtools(args.samtool)

    genome = library.filter("hs37d5", Source.EBI_ALT)

    downloader = Downloader()
    if downloader.need_download(genome):
        downloader.download(genome)
    
    decompressor = Decompressor(samtools)
    if decompressor.need_decompression(genome):
        decompressor.decompress(genome)

    compressor = Compressor(samtools)
    if compressor.need_compression(genome):
        compressor.compress(genome)

    t = samtools.get_file_type(genome.initial_name)
    t = samtools.fasta_index(genome.initial_name)
    t = samtools.make_dictionary(genome.initial_name)