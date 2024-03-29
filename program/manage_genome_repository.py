import argparse
import logging
from pathlib import Path
from reference.genome_repository import GenomeRepository
import pprint


def download(args):
    repository: GenomeRepository = GenomeRepository.build(
        args.csv, args.root, args.external
    )
    genomes = repository.filter(args.id)
    logging.debug(f"Found {len(genomes)} genomes matching these criterias:")
    logging.debug(pprint.pformat(genomes))
    unique_genomes = dict()
    if not args.all:
        for genome in genomes:
            if genome.code not in unique_genomes:
                unique_genomes[genome.code] = genome

    if len(genomes) != len(unique_genomes):
        delta = len(genomes) - len(unique_genomes)
        logging.info(
            f"Found {len(genomes)} genomes matching the criteria(s) but filtered out {delta} "
            "because they were duplicate of the same genome coming from different sources. "
            "To override this behavior specify --all."
        )
    
    logging.info("Genomes to download:")
    for index, genome in enumerate(unique_genomes.values()):
        logging.info(f"#{index} {str(genome)}")

    for genome in unique_genomes.values():
        logging.info(f"Adding to library: {str(genome)}")
        repository.add(genome)

def list_(args):
    repository: GenomeRepository = GenomeRepository.build(
        args.csv, args.root, args.external
    )
    


def add_(args):
    repository: GenomeRepository = GenomeRepository.build(
        args.csv, args.root, args.external
    )


def delete_(args):
    repository: GenomeRepository = GenomeRepository.build(
        args.csv, args.root, args.external
    )
    genomes = repository.filter(args.id)
    if not args.yes:
        y = "n"
        while y.lower() != "y":
            try:
                y = input(f"You're about to delete {len(genomes)} from repository. Press Y to continue.")
            except EOFError:
                logging.info("Aborting.")
                return
    for genome in genomes:
        repository.delete(genome)



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
        "--external",
        type=str,
        metavar="path",
        help="Path to the folder containing external executables.",
        default=None,
    )
    parser.add_argument(
        "--verbose",
        help="Set logging level",
        choices=["INFO", "DEBUG", "WARNING", "ERROR"],
        type=str,
        default="INFO",
    )

    action_subparser = parser.add_subparsers(title="subcommands")
    download_action = action_subparser.add_parser(
        "download",
        help="Download a reference genome file",
        aliases=["dl"],
    )
    download_action.set_defaults(func=download)
    download_action.add_argument(
        "--id", type=str, help="ID of the reference genome to download."
    )
    download_action.add_argument(
        "--all",
        action="store_true",
        help="Download every genome matching the criteria, even if it's a duplicate.",
    )

    list_action = action_subparser.add_parser(
        "list", help="List available reference genomes.", aliases=["ls"]
    )
    list_action.add_argument(
        "--on-disk", type=str, help="ID of the reference genome to download."
    )
    list_action.set_defaults(func=list_)

    add_action = action_subparser.add_parser(
        "add", help="Add a reference genome to the repository."
    )
    add_action.set_defaults(func=add_)
    delete_action = action_subparser.add_parser(
        "delete", help="Delete a reference genome from disk.", aliases=["rm"]
    )   
    delete_action.add_argument(
        "--yes", action="store_true", help="Skip confirmation."
    )
    delete_action.set_defaults(func=delete_)

    #args = parser.parse_args(["download", "--id", "hs37d5"])
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.__dict__[args.verbose])
    
    if "func" in args.__dict__:
        args.func(args)
    else:
        parser.print_help()
