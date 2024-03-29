#!/usr/bin/env python3
# coding: utf8
#
# Part of the WGS Extract (https://wgse.bio/) system (standalone)
#
# Copyright (C) 2022-2024 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""
Process compressed FASTA files to generate statistics on the presence of unknown bases (literal 'N's inside the FASTA).

Standalone script to process a reference model FASTA to determine the inclusion of N's in the base pair sequence.
Reads the compressed FASTA. Takes as a parameter the FASTA file. Also expects and reads in the DICT file determined
from the FASTA file name -- easier to do one pass if we have the sequence lengths.  Processes all sequences.
Detailed N segment stats writtem to RefModel_nbin.csv. Key stats per sequence and totals written to RefModel_nbuc.csv.
Also creates nat log bucket counts of number of N's; 1000 buckets per sequence. Per David Vance's display criteria.
Finally, a variation of the _nbin.csv file is created in BED format as _nreg.bed.

Sequences in a reference model are arbitrarily broken up into short lines of 50-80 characters (base-pairs) each. All
white space is in-material and to be ignored.  We want to not just count how many N values in a sequence but also
(1) how many regions of contigious N's (with their mean and std dev of size) and
(2) a thousand buckets for each sequence representing the round(natlog(# of N's)) in that bucket.
Sequence of contigous N values can overlap the short line and bucket boundaries and must be handled appropriately.
Due to rounding errors, the last bucket may not be the same size of base-pairs as the rest of the buckets.

Detects if a FASTQ instead of FASTA and reports error. Skips any sequence where the name is not found in the
DICT file read in to start.  Key is not just the SN but also LN fields in the DICT file.
"""

import argparse
import math
import pathlib
import logging
import time
import typing
import statistics

from .buckets import Buckets
from .fasta_file import FastaFile, Sequence
from .external import External
from .genome import Genome


class NStatisticsFiles:
    """Manage the creation of statistics files for Ns."""

    def __init__(
        self, 
        fasta_file: FastaFile,
        long_run_threshold: int = 300,
        buckets_number: int = 1000
    ) -> None:
        self._long_run_threshold = long_run_threshold
        self._buckets_number = buckets_number
        self._fasta_file = fasta_file

    def get_nbin(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract runs of N: BIN definition file\n",
            f"#Processing Ref Model: {self._fasta_file.model_name} with >{self._long_run_threshold}bp of N runs\n",
            "#SN\tBinID\tStart\tSize\n",
        ]
        for sequence in sequences:
            for index, run in enumerate(
                sequence.filter(lambda x: x.length > self._long_run_threshold)
            ):
                row = f"{sequence.name}\t{index+1}\t{run.start:,}\t{run.length:,}\n"
                lines.append(row)
        return lines

    def get_bed(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract runs of N: BED file of bin definitions\n",
            f"#Processing Ref Model: {self._fasta_file.model_name} with >{self._long_run_threshold}bp of N runs\n",
            "#SN\tStart\tStop\n",
        ]

        for sequence in sequences:
            for run in sequence.filter(lambda x: x.length > self._long_run_threshold):
                row = f"{sequence.name}\t{run.start}\t{run.start+run.length}\n"
                lines.append(row)
        return lines

    def get_nbuc(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract generated Sequence of N summary file\n",
            f"#Model {self._fasta_file.model_name} with >{self._long_run_threshold}bp Sequence of N and {self._buckets_number} buckets per sequence\n",
            "#Seq\tNumBP\tNumNs\tNumNreg\tNregSizeMean\tNregSizeStdDev\tSmlNreg\tBuckSize\tBucket Sparse List (bp start, ln value) when nonzero\n",
        ]
        for sequence in sequences:
            bucket_lenght = int(math.floor(sequence.length / self._buckets_number))

            long_runs = sequence.filter(lambda x: x.length > self._long_run_threshold)
            short_runs = sequence.filter(lambda x: x.length <= self._long_run_threshold)

            long_run_lengths = [x.length for x in long_runs]
            short_run_lenghts = [x.length for x in short_runs]

            long_runs_total = sum(long_run_lengths)
            short_runs_total = sum(short_run_lenghts)

            long_run_avg = 0
            long_run_stdev = 0

            if len(long_run_lengths) > 1:
                long_run_avg = int(statistics.mean(long_run_lengths))
                long_run_stdev = int(statistics.stdev(long_run_lengths))

            row = f"{sequence.name}\t{sequence.length}\t{long_runs_total}\t{len(long_run_lengths)}\t{long_run_avg}\t{long_run_stdev}\t{short_runs_total}\t{bucket_lenght}"

            buckets = Buckets(
                sequence, self._buckets_number, self._long_run_threshold
            )
            for index, bucket_count in buckets.buckets.items():
                val = round(math.log(bucket_count) if bucket_count > 1 else 0)
                start = index * bucket_lenght

                if val > 0:
                    row += f"\t{start}\t{val}"
            row += "\n"
            lines.append(row)
        return lines

    def _generate_file(self, path: pathlib.Path, lines: typing.List[str]):
        with path.open("wt", encoding="utf8") as f:
            f.writelines(lines)

    def _generate_files(self, sequences: typing.Dict[str, typing.List[Sequence]]):
        self._generate_file(self._fasta_file.genome.nbuc, self.get_nbuc(sequences))
        self._generate_file(self._fasta_file.genome.bed, self.get_bed(sequences))
        self._generate_file(self._fasta_file.genome.nbin, self.get_nbin(sequences))

    def generate_stats(self):
        logging.info(f"{self._fasta_file.genome.final_name.name}: Counting Ns.")
        sequences = self._fasta_file.count_letters("N")
        logging.info(f"{self._fasta_file.genome.final_name.name}: Finished counting Ns.")
        self._generate_files(sequences)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("unknown_bases_stats.py")
    parser.add_argument(
        "--reference",
        help="Indicate the path of the bgzip compressed reference genome",
        type=pathlib.Path,
        default=pathlib.Path("reference\\genomes\\_hs37d5.fa\\seq1.fa.gz"),
    )
    parser.add_argument(
        "--reference-dict",
        help="Indicate the path of dictionary for the reference genome",
        type=pathlib.Path,
        default=None,
    )
    parser.add_argument(
        "--long-run-threshold",
        help="Indicate the treshold for long runs to be considered when making .buc files.",
        type=int,
        default=300,
    )
    parser.add_argument(
        "--buckets-number",
        help="Indicate the total number of buckets to be considered when making .buc files.",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--external",
        help="Indicate the root directory for 3rd parties executables (only needed to generate the .dict file if not already available).",
        type=pathlib.Path,
        default=None,
    )
    parser.add_argument("--verbose", help="Set logging level", choices=["INFO", "DEBUG", "WARNING", "ERROR"], type=str, default="INFO")
    args = parser.parse_args()
    
    logging.getLogger().setLevel(logging.__dict__[args.verbose])

    start = time.time()
    genome = Genome(args.reference)
    
    if not genome.dict.exists():
        external = External(args.external)
        external.make_dictionary(args.reference, genome.dict)
    
    fasta_file = FastaFile(genome)
    unknown_bases_stats = NStatisticsFiles(
        fasta_file, args.long_run_threshold, args.buckets_number,
    )
    unknown_bases_stats.generate_stats()
    end = time.time()
    logging.info(f"Time {end-start}s")
