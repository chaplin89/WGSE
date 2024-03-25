import logging
from pathlib import Path
import sys
from reference.decompressor import Decompressor
from reference.compressor import Compressor
from reference.downloader import Downloader
from reference.file_type_checker import FileTypeChecker
from reference.genome_repository import GenomeRepository
from reference.samtools import Samtools
from reference.csv_metadata_loader import CsvMetadataLoader
import argparse

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    
    if "win" in sys.platform:
        default_samtools = Path("cygwin64", "usr", "local", "bin")

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
        default=default_samtools,
    )
    
    parser.add_argument(
        "--verbosity",
        type=int,
        metavar="integer",
        help="Set the verbosity level.",
        default=10,
    )

    action_subparser = parser.add_subparsers(title="subcommands")
    download_action = action_subparser.add_parser(
        "download", help="Download a reference genome file", aliases=["dl"]
    )
    list_action = action_subparser.add_parser(
        "list", help="List available reference genomes.", aliases=["ls"]
    )
    add_action = action_subparser.add_parser(
        "add", help="Add a reference genome to the repository."
    )
    delete_action = action_subparser.add_parser(
        "delete", help="List available reference genomes.", aliases=["rm"]
    )
    
    download_action.add_argument(
        "id", type=str, help="ID of the reference genome to download."
    )
    list_action.add_argument(
        "--on-disk", type=str, help="ID of the reference genome to download."
    )

    args = parser.parse_args(["download","hg37"])
    parser.print_help()
    
    samtools = Samtools(args.samtools)
    file_type_checker = FileTypeChecker(samtools)
    downloader = Downloader()
    decompressor = Decompressor(file_type_checker)
    compressor = Compressor(samtools)
    metadata = CsvMetadataLoader(args.csv)

    repository = GenomeRepository(
        metadata.genomes,
        args.root,
        file_type_checker,
        downloader,
        compressor,
        decompressor,
    )

    for genome in metadata.genomes:
        repository.add_to_library(genome)

    # genome = library.filter("hs37d5", Source.EBI_ALT)
