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
from math import floor
import math
from pathlib import Path
import re
import gzip
import logging
import time
from typing import Dict, List, TextIO
from collections import OrderedDict
from statistics import mean, stdev


class Run:
    """Represent a single run of Ns"""

    def __init__(self) -> None:
        self.start = None
        self.length = None

    def open(self, position: int):
        """_summary_

        Args:
            position (int): _description_
        """
        self.start = position

    def close(self, length: int):
        """_summary_

        Args:
            length (int): _description_
        """
        self.length = length


class Sequence:
    """Represent a sequence containing multiple runs of N."""

    def __init__(self, name: str, length: int) -> None:
        """_summary_

        Args:
            name (str): _description_
            length (int): _description_
        """
        self.runs: List[Run] = []
        self.name: str = name
        self.length = length
        self._current_run: Run = None

    def is_run_open(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        return self._current_run is not None

    def open_run(self, position: int) -> None:
        """Open a new run of Ns in this sequence.

        Args:
            position (int): Position where the run is starting.

        Raises:
            RuntimeError: Opening a sequence that is already opened.
        """
        if self._current_run is not None:
            raise RuntimeError("Trying to open a an already opened run of Ns.")
        # Ensure no overlapping runs
        if len(self.runs) > 0:
            previous_end = self.runs[-1].start + self.runs[-1].length
            if position < previous_end + 1:
                raise RuntimeError("Trying to open a run overlapping with previous runs.")
        
        self._current_run = Run()
        self._current_run.open(position)

    def filter_runs(self, length_greater_than: int):
        return [x for x in self.runs if x.length > length_greater_than]
    
    def split_in_buckets(self, buckets_number : int) -> OrderedDict[int, int]:
        if buckets_number > self.length:
            raise RuntimeError(f"Number of buckets should be at least {self.length}")
        
        buckets = OrderedDict()
        bucket_size = self.length / buckets_number
        
        for run in self.runs:
            start = run.start
            end = run.start + run.length
            index_bucket_start =  floor(start / bucket_size)
            index_bucket_end = floor(end / bucket_size)
            
            for bucket in range(index_bucket_start, index_bucket_end + 1):
                start_bucket_offset = max(bucket * bucket_size, start)
                end_bucket_offset = min(start_bucket_offset + bucket_size - 1, end)
                if bucket not in buckets:
                    buckets[bucket] = 0
                runs_count = end_bucket_offset - start_bucket_offset
                buckets[bucket] += runs_count + 1
        return buckets

    def close_run(self, position: int) -> None:
        """Close an open run of Ns in this sequence.

        Args:
            position (int): Position where the run is ending.

        Raises:
            RuntimeError: Closing a sequence that is not open or closing a
            sequence with a position smaller than the open.
        """
        if self._current_run is None:
            raise RuntimeError("Trying to close an already closed run of Ns.")
        if position <= self._current_run.start:
            raise RuntimeError(
                f"Trying to close a run of Ns with position {position} and start {self._current_run.start}"
            )
        run_length = position - self._current_run.start
        self._current_run.close(run_length)
        self.runs.append(self._current_run)
        self._current_run = None


class Stats:
    def __init__(self, long_threshold: int, buckets_number, model: str) -> None:
        self._long_threshold = long_threshold
        self._buckets_number = buckets_number
        self._model = model

    def get_nbin(self, sequences: List[Sequence]):
        lines = [
            f"#WGS Extract runs of N: BIN definition file\n",
            f"#Processing Ref Model: {self._model} with >{self._long_threshold}bp of N runs\n",
            f"#SN\tBinID\tStart\tSize\n",
        ]
        for sequence in sequences:
            for index, run in enumerate(sequence.filter_runs(self._long_threshold)):
                row = f"{sequence.name}\t{index+1}\t{run.start:,}\t{run.length:,}\n"
                lines.append(row)
        return lines

    def get_bed(self, sequences: List[Sequence]):
        lines = [
            f"#WGS Extract runs of N: BED file of bin definitions\n",
            f"#Processing Ref Model: {self._model} with >{self._long_threshold}bp of N runs\n",
            f"#SN\tStart\tStop\n",
        ]

        for sequence in sequences:
            for run in sequence.filter_runs(self._long_threshold):
                row = f"{sequence.name}\t{run.start}\t{run.start+run.length}\n"
                lines.append(row)
        return lines

    def get_nbuc(self, sequences: List[Sequence]):
        lines = [
            f"#WGS Extract generated Sequence of N summary file\n",
            f"#Model {self._model} with >{self._long_threshold}bp Sequence of N and {self._buckets_number} buckets per sequence\n",
            f"#Seq\tNumBP\tNumNs\tNumNreg\tNregSizeMean\tNregSizeStdDev\tSmlNreg\tBuckSize\tBucket Sparse List (bp start, ln value) when nonzero\n",
        ]
        for sequence in sequences:
            long_runs = sequence.filter_runs(self._long_threshold)
            lengths = [x.length for x in long_runs]
            average_length = mean(lengths) if len(lengths) > 0 else 0
            standard_deviation = stdev(lengths) if len(lengths) > 1 else 0
            total = sum(lengths)
            
            row = f"{sequence.name}\tNumBP\t{total}\t{len(lengths)}\t{average_length}\t{standard_deviation}\tSmlNreg\tBuckSize\t\n"
            
            for index, bucket in range(sequence.split_in_buckets(self._buckets_number)):
                val = round(math.log(bucket) if bucket > 1 else 0)
                # changed to print as sparse list as most of the 1000 entries are zero; so now bucket start (bp) then len (Ncnt)
                if val > 0:
                    start = index * bucksize
                    lines.append(f'\t{start}\t{val}')
            
            lines.append(row)
        return lines

class UnknownBasesStats:
    def __init__(
        self, path: Path, long_run_threshold: int = 300, buckets_number: int = 1000
    ):
        self._long_run_threshold = long_run_threshold
        self._buckets_number = buckets_number
        self._path = path

    def _process_file(self, fp: TextIO) -> List[Sequence]:
        sequences = OrderedDict()
        pattern = re.compile(r"N+")

        position = None
        current_sequence = None

        for line in fp:
            line = line.rstrip("\n")

            # Comment: skip
            if len(line) == 0 or line[0] == "#":
                continue

            # Start of a new sequence
            if line[0] == "+":
                raise RuntimeError(
                    f"Expected a FASTA reference model, got a FASTQ sequencer data."
                )

            if line[0] == ">":
                if current_sequence is not None and current_sequence.is_run_open():
                    current_sequence.close_run(position)
                sequence_name = line.split()[0][1:]
                sequence_length = int(line.split()[2].split(':')[4])
                if len(sequences) > 0:
                    pass
                    # return list(sequences.values())
                logging.info(f"Processing sequence {sequence_name}")
                if sequence_name in sequences:
                    raise RuntimeError(f"Found a duplicated sequence: {sequence_name}")
                current_sequence = Sequence(sequence_name, sequence_length)
                sequences[sequence_name] = current_sequence
                position = 0
                continue

            matches = list(pattern.finditer(line))

            # Next run found but there's still an open run: close
            if len(matches) == 0:
                if current_sequence.is_run_open():
                    current_sequence.close_run(position)

            for match in matches:
                if match.start() == 0:
                    # If open, keep open: is continuing from the previous line.
                    if not current_sequence.is_run_open():
                        current_sequence.open_run(position)
                else:
                    if current_sequence.is_run_open():
                        current_sequence.close_run(position)
                    current_sequence.open_run(match.start() + position)

                if match.end() < len(line):
                    current_sequence.close_run(match.end() + position)
            position += len(line)

        # File was terminated but there's still an open run: close
        if current_sequence.is_run_open():
            current_sequence.close_run(position)
        return list(sequences.values())

    def _count_unknown_bases(self):
        logging.info(f"Starting to process file: {self._path.name}")
        with gzip.open(self._path, "rt") as f:
            return self._process_file(f)

    def _generate_file(self, suffix: str, lines: List[Sequence]):
        target_name = self._path.stem.strip("".join(self._path.suffixes)) + suffix
        target = self._path.parent.joinpath(target_name)
        with open(target, "wt") as f:
            f.writelines(lines)

    def _generate_files(self, unknown_bases: Dict[str, List[Sequence]]):
        stats = Stats(self._long_run_threshold, self._buckets_number, self._path.name)
        self._generate_file("_nbuc.csv", stats.get_nbuc(unknown_bases))
        self._generate_file(".bed", stats.get_bed(unknown_bases))
        self._generate_file("_nbin.csv", stats.get_nbin(unknown_bases))

    def get_stats(self):
        unknown_bases = self._count_unknown_bases()
        self._generate_files(unknown_bases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("unknown_bases_stats.py")
    parser.add_argument(
        "--reference",
        help="Indicate the path of the bgzip compressed reference genome",
        type=Path,
        default=Path("reference\\genomes\\hs37d5.fa.gz"),
    )
    parser.add_argument(
        "--verbose", help="Set logging level [1,5]", type=int, default=2
    )
    args = parser.parse_args()
    logging.getLogger().setLevel(args.verbose)

    start = time.time()
    unknown_bases_stats = UnknownBasesStats(args.reference)
    unknown_bases_stats.get_stats()
    end = time.time()
    logging.info(f"Time {end-start}s")
