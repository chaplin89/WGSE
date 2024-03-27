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
import sys
import time
import typing
import collections
import statistics
from .fasta_file import FastaFile, Sequence
from .samtools import Samtools
from .genome import Genome


class Buckets:
    """Represent a sequence partitioned into subsequences of fixed lenght called buckets"""

    def __init__(
        self, sequence: Sequence, buckets_number: int, long_run_threshold: int
    ) -> None:
        self._buckets_number = buckets_number
        self._long_run_threshold: int = long_run_threshold
        self._sequence: Sequence = sequence
        self.buckets = self._make_buckets()

    def _make_buckets(self) -> typing.OrderedDict[int, int]:
        if self._buckets_number > self._sequence.length:
            # Can't have less than 1 N per bucket.
            return collections.OrderedDict()
            raise ValueError(
                f"Number of buckets should be at least {self._sequence.length}"
            )

        buckets = collections.OrderedDict()
        bucket_size = int(math.floor(self._sequence.length / self._buckets_number))

        for run in self._sequence.filter(
            lambda x: x.length >= self._long_run_threshold
        ):
            # Determine how many buckets this run is spanning
            start = run.start
            end = run.start + run.length
            index_bucket_start = int(math.floor(start / bucket_size))
            index_bucket_end = int(math.floor(end / bucket_size))

            # Iterate over each bucket determining how many
            for bucket in range(index_bucket_start, index_bucket_end + 1):
                start_bucket_offset = max(bucket * bucket_size, start)
                end_bucket_offset = min(bucket * bucket_size + bucket_size - 1, end)
                runs_count = end_bucket_offset - start_bucket_offset
                if runs_count < 1:
                    continue
                if bucket not in buckets:
                    buckets[bucket] = 0
                buckets[bucket] += runs_count
        return buckets

    def __getitem__(self, key):
        return self.buckets[key]


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
            ).buckets
            for index, bucket_count in buckets.items():
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

    def _generate_files(self, unknown_bases: typing.Dict[str, typing.List[Sequence]]):
        self._generate_file(self._fasta_file.genome.nbuc, self.get_nbuc(unknown_bases))
        self._generate_file(self._fasta_file.genome.bed, self.get_bed(unknown_bases))
        self._generate_file(self._fasta_file.genome.nbin, self.get_nbin(unknown_bases))

    def generate_stats(self):
        sequences = self._fasta_file.split_into_sequences()
        self._generate_files(sequences)


if __name__ == "__main__":
    if "win" in sys.platform:
        default_samtools = pathlib.Path("cygwin64", "usr", "local", "bin")
    else: 
        default_samtools = None
        
    
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
        "--samtools",
        help="Indicate the root directory for samtools (only needed to generate the .dict file if not already available).",
        type=pathlib.Path,
        default=default_samtools,
    )
    parser.add_argument("--verbose", help="Set logging level [1-5]", type=int, default=1)
    args = parser.parse_args()
    
    logging.getLogger().setLevel(args.verbose*10)

    start = time.time()
    genome = Genome(args.reference)
    if not genome.dict.exists():
        samtools = Samtools(args.samtools)
        samtools.make_dictionary(args.reference, genome.dict)
    
    unknown_bases_stats = NStatisticsFiles(
        FastaFile(genome), args.long_run_threshold, args.buckets_number,
    )
    unknown_bases_stats.generate_stats()
    end = time.time()
    logging.info(f"Time {end-start}s")
