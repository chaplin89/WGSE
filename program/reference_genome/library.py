import csv
import enum
from pathlib import Path
from typing import List, Union, Tuple
from decompressor import Decompressor
from compressor import Compressor
from downloader import Downloader
from file_type_checker import FileTypeChecker, Type
from genome import Genome, Source
from samtools import Samtools
from csv_metadata_loader import CsvMetadataLoader
import argparse


class GenomeRepository:

    def __init__(
        self,
        reference_genomes: List[Genome],
        library_path: Path,
        type_checker: FileTypeChecker,
        downloader: Downloader,
        compressor: Compressor,
        decompressor: Decompressor,
    ) -> None:
        """Initialize the object

        Args:
            reference_genomes (List[Genome]): List of known reference genomes
            library_path (Path): Root directory for the repository
            downloader (Downloader): An instance of Downloader
            compressor (Compressor): An instance of Compressor
            decompressor (Decompressor): An instance of Decompressor
        """
        self._reference_genomes = reference_genomes
        self._library_path = library_path
        self._type_checker = type_checker
        self._downloader = downloader
        self._decompressor = decompressor
        self._compressor = compressor
        self._rebase_reference_genomes(self._reference_genomes, self._library_path)

    def determine_target_name(self, file: Path):
        if file.suffix in [".bz2", ".bz", ".zip", ".7z"]:
            gz_compressed = Path(file.parent, file.stem + ".gz")
        elif file.suffix != ".gz":
            gz_compressed = Path(str(file) + ".gz")
        else:
            gz_compressed = file

        return gz_compressed

    def _rebase_reference_genomes(self, genomes: List[Genome], root: Path):
        """Fix the filename for reference genomes so they're pointing to a file
        under the root directory of this repository.

        Args:
            genomes (List[Genome]): List of reference genomes
            root (Path): Root path for the library
        """
        for genome in genomes:
            genome.file = root.joinpath(genome.file)

    def to_bgzip(self, genome: Genome):
        type = self._type_checker.get_type(genome.file)

        if type == Type.BGZIP:
            return
        
        if type != Type.DECOMPRESSED:
            decompressed = self._decompressor.decompress(genome)
            genome.file.unlink()
            genome.file = decompressed
        
        compressed = self._compressor.compress(genome)
        genome.file.unlink()
        genome.file = compressed

    def add_to_library(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        if self._downloader.need_download(genome):
            self._downloader.download(genome, None)

        target = self.determine_target_name(genome.file)
        if target.exists():
            genome.file = target
        
        type = self._type_checker.get_type(genome.file)

        if type != Type.BGZIP:
            try:
                self.to_bgzip(genome)
            except:
                pass

        #t = samtools.fasta_index(genome.file)
        #t = samtools.make_dictionary(genome.file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Manage a library of reference genomes.",
        prog=Path(__file__).name,
        usage="%(prog)s [options]",
    )

    parser.add_argument(
        "--csv",
        type=Path,
        metavar="path",
        help="A CSV file containing a description on known reference genomes.",
        default=Path("reference", "seed_genomes.csv"),
    )

    parser.add_argument(
        "--root",
        type=Path,
        metavar="path",
        help="Path to the root folder of the library.",
        default=Path("reference", "genomes"),
    )

    parser.add_argument(
        "--samtools",
        type=str,
        metavar="path",
        help="Path to the folder containing samtool executables.",
        default=Path("cygwin64", "usr", "local", "bin"),
    )

    download_subparser = parser.add_subparsers(help="Download a reference genome file")
    download_action = download_subparser.add_parser(
        "download", help="Download a reference genome file"
    )
    download_action.add_argument(
        "id", type=str, help="ID of the reference genome to download."
    )

    args = parser.parse_args()
    parser.print_help()

    samtools = Samtools(args.samtools)
    file_type_checker = FileTypeChecker(samtools)
    downloader = Downloader()
    decompressor = Decompressor(file_type_checker)
    compressor = Compressor(file_type_checker)
    metadata = CsvMetadataLoader(args.csv)

    repository = GenomeRepository(
        metadata.genomes, args.root, file_type_checker, downloader, compressor, decompressor
    )

    for genome in metadata.genomes:
        repository.add_to_library(genome)

    # genome = library.filter("hs37d5", Source.EBI_ALT)