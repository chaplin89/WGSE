import csv
import enum
import logging
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
    """Manage a repository of reference genome files"""

    def __init__(
        self,
        reference_genomes: List[Genome],
        library_path: Path,
        type_checker: FileTypeChecker,
        downloader: Downloader,
        compressor: Compressor,
        decompressor: Decompressor,
    ) -> None:
        self._reference_genomes = reference_genomes
        self._library_path = library_path
        self._type_checker = type_checker
        self._downloader = downloader
        self._decompressor = decompressor
        self._compressor = compressor
        self._rebase_reference_genomes(self._reference_genomes, self._library_path)

    def determine_target_name(self, file: Path):
        if file.suffix in [".bz2", ".bz", ".zip", ".7z"]:
            gz_compressed = file.with_suffix(".gz")
        elif file.suffix != ".gz":
            gz_compressed = Path(str(file) + ".gz")
        else:
            gz_compressed = file

        return gz_compressed

    def _rebase_reference_genomes(self, genomes: List[Genome], root: Path):
        for genome in genomes:
            genome.initial_name = root.joinpath(genome.initial_name)

    def to_bgzip(self, genome: Genome):
        type = self._type_checker.get_type(genome.initial_name)

        if type == Type.RAZF_GZIP:
            logging.info(f"Determined type RAZF for {genome.initial_name}. Skipping.")
            # TODO: figure out how to decompress this format.
            return

        if type == Type.BGZIP:
            logging.info(f"Determined type BGZIP for {genome.initial_name}. Skipping.")
            return
       
        if type != Type.DECOMPRESSED:            
            logging.info(f"Determined type {type} for {genome.initial_name}. Decompressing.")
            decompressed = self._decompressor.decompress(genome)
            genome.initial_name.unlink()
            genome.initial_name = decompressed
        
        logging.info(f"Compressing {genome.initial_name} into bgzip.")
        compressed = self._compressor.compress(genome)
        genome.initial_name.unlink()
        genome.initial_name = compressed

    def add_to_library(self, genome: Genome):
        """Add a genome to the library

        Args:
            genome (Genome): _description_
        """
        
        logging.info(f"Adding Genome {genome.code} to library.")
        if self._downloader.need_download(genome):
            logging.info(f"Downloading Genome {genome.code}.")
            self._downloader.download(genome, None)

        target = self.determine_target_name(genome.initial_name)
        logging.info(f"Genome target name: {target}.")

        if target.exists():
            genome.initial_name = target
        
        type = self._type_checker.get_type(genome.initial_name)
        logging.info(f"Genome file type determined: {type}.")

        if type != Type.BGZIP:
            try:
                logging.info(f"Starting conversion of {type} to bgzip.")
                self.to_bgzip(genome)
            except Exception as e:
                logging.error(e)

        #t = samtools.fasta_index(genome.file)
        #t = samtools.make_dictionary(genome.file)


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

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
    compressor = Compressor(samtools)
    metadata = CsvMetadataLoader(args.csv)

    repository = GenomeRepository(
        metadata.genomes, args.root, file_type_checker, downloader, compressor, decompressor
    )

    for genome in metadata.genomes:
        repository.add_to_library(genome)

    # genome = library.filter("hs37d5", Source.EBI_ALT)