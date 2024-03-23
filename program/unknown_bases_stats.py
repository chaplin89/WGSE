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
import re
import gzip
import logging
import time
import typing
import collections
import statistics


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
        if length <= 0:
            raise IndexError("Length should be greater than zero.")
        self.runs: typing.List[Run] = []
        self.name: str = name
        self.length: int = length
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
            RuntimeError: Opening a sequence that is already opened,
            opening a run with a negative position,
            opening a run that is overlapping with the previous.
        """
        if self._current_run is not None:
            raise RuntimeError("Trying to open a an already opened run of Ns.")
        if position < 0:
            raise IndexError("Trying to open a run with a negative position.")
        # Ensure no overlapping runs
        if len(self.runs) > 0:
            previous_end = self.runs[-1].start + self.runs[-1].length
            if position < previous_end + 1:
                raise ValueError("Trying to open a run overlapping with previous runs.")

        self._current_run = Run()
        self._current_run.open(position)

    def filter(self, criteria: typing.Callable[[Run], bool]):
        return [x for x in self.runs if criteria(x)]

    def split_in_buckets(self, buckets_number: int) -> typing.OrderedDict[int, int]:
        """_summary_

        Args:
            buckets_number (int): _description_

        Raises:
            RuntimeError: _description_

        Returns:
            OrderedDict[int, int]: _description_
        """
        if buckets_number > self.length:
            raise ValueError(f"Number of buckets should be at least {self.length}")

        buckets = collections.OrderedDict()
        bucket_size = int(math.floor(self.length / buckets_number))

        for run in self.runs:
            start = run.start
            end = run.start + run.length
            index_bucket_start = int(math.floor(start / bucket_size))
            index_bucket_end = int(math.floor(end / bucket_size))

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
        if position > self.length:
            raise RuntimeError(
                "Trying to close a run with a position that is greater than the whole length of the sequence."
            )
        run_length = position - self._current_run.start
        self._current_run.close(run_length)
        self.runs.append(self._current_run)
        self._current_run = None

    def end(self, position: int) -> None:
        pass


class Stats:
    def __init__(self, long_threshold: int, buckets_number, model: str) -> None:
        self._long_threshold = long_threshold
        self._buckets_number = buckets_number
        self._model = model

    def get_nbin(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract runs of N: BIN definition file\n",
            f"#Processing Ref Model: {self._model} with >{self._long_threshold}bp of N runs\n",
            "#SN\tBinID\tStart\tSize\n",
        ]
        for sequence in sequences:
            for index, run in enumerate(
                sequence.filter(lambda x: x.length > self._long_threshold)
            ):
                row = f"{sequence.name}\t{index+1}\t{run.start:,}\t{run.length:,}\n"
                lines.append(row)
        return lines

    def get_bed(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract runs of N: BED file of bin definitions\n",
            f"#Processing Ref Model: {self._model} with >{self._long_threshold}bp of N runs\n",
            "#SN\tStart\tStop\n",
        ]

        for sequence in sequences:
            for run in sequence.filter(lambda x: x.length > self._long_threshold):
                row = f"{sequence.name}\t{run.start}\t{run.start+run.length}\n"
                lines.append(row)
        return lines

    def get_nbuc(self, sequences: typing.List[Sequence]):
        lines = [
            "#WGS Extract generated Sequence of N summary file\n",
            f"#Model {self._model} with >{self._long_threshold}bp Sequence of N and {self._buckets_number} buckets per sequence\n",
            "#Seq\tNumBP\tNumNs\tNumNreg\tNregSizeMean\tNregSizeStdDev\tSmlNreg\tBuckSize\tBucket Sparse List (bp start, ln value) when nonzero\n",
        ]
        for sequence in sequences:
            bucket_lenght = int(math.floor(sequence.length / self._buckets_number))

            long_runs = sequence.filter(lambda x: x.length > self._long_threshold)
            short_runs = sequence.filter(lambda x: x.length <= self._long_threshold)

            long_run_lengths = [x.length for x in long_runs]
            short_run_lenghts = [x.length for x in short_runs]

            long_run_avg = 0
            long_run_stdev = 0

            if len(long_run_lengths) > 1:
                long_run_avg = int(statistics.mean(long_run_lengths))
                long_run_stdev = int(statistics.stdev(long_run_lengths))

            long_runs_total = sum(long_run_lengths)
            short_runs_total = sum(short_run_lenghts)

            row = f"{sequence.name}\t{sequence.length}\t{long_runs_total}\t{len(long_run_lengths)}\t{long_run_avg}\t{long_run_stdev}\t{short_runs_total}\t{bucket_lenght}"

            buckets = sequence.split_in_buckets(self._buckets_number).items()
            for index, bucket_count in buckets:
                val = round(math.log(bucket_count) if bucket_count > 1 else 0)
                start = index * bucket_lenght

                if val > 0:
                    row += f"\t{start}\t{val}"
            row += "\n"
            lines.append(row)
        return lines


class UnknownBasesStats:
    def __init__(
        self,
        path: pathlib.Path,
        long_run_threshold: int = 300,
        buckets_number: int = 1000,
    ):
        self._long_run_threshold = long_run_threshold
        self._buckets_number = buckets_number
        self._path = path

    def _process_file(self, fp: typing.TextIO) -> typing.List[Sequence]:
        sequences = collections.OrderedDict()
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
                    "Expected a FASTA reference model, got a FASTQ sequencer data."
                )

            if line[0] == ">":
                if current_sequence is not None and current_sequence.is_run_open():
                    current_sequence.close_run(position)
                sequence_name = line.split()[0][1:]
                sequence_length = int(line.split()[2].split(":")[4])
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
        current_sequence.end(position)
        return list(sequences.values())

    def _count_unknown_bases(self):
        logging.info(f"Starting to process file: {self._path.name}")
        with gzip.open(self._path, "rt") as f:
            return self._process_file(f)

    def _generate_file(self, suffix: str, lines: typing.List[Sequence]):
        target_name = self._path.stem.strip("".join(self._path.suffixes)) + suffix
        target = self._path.parent.joinpath(target_name)
        with open(target, "wt", encoding="utf8") as f:
            f.writelines(lines)

    def _generate_files(self, unknown_bases: typing.Dict[str, typing.List[Sequence]]):
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
        type=pathlib.Path,
        default=pathlib.Path("reference\\genomes\\seq1.txt.gz"),
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
