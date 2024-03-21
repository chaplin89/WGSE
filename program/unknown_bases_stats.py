#!/usr/bin/env python3
# coding: utf8
#
# Counting Reference Model Final Assembly N's (BED, region, etc output files)
#
# Part of the WGS Extract (https://wgse.bio/) system (standalone)
#
# Copyright (C) 2022-2024 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""
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
from pathlib import Path
import sys
import os
import re
import math
import gzip
import logging
import enum
import time
from typing import Dict, List, TextIO, Tuple


class Sequence:
    """Represent a sequence containing multiple runs of N."""

    def __init__(self, name: str) -> None:
        self.runs = []
        self._current_start = None
        self._name = name

    def is_open(self) -> bool:
        return self._current_start is not None

    def open(self, position) -> None:
        """_summary_

        Args:
            position (_type_): _description_

        Raises:
            RuntimeError: _description_
        """
        if self._current_start is not None:
            raise RuntimeError("Trying to open a an already opened run of Ns.")
        self._current_start = position

    def close(self, position) -> None:
        """_summary_

        Args:
            position (_type_): _description_

        Raises:
            RuntimeError: _description_
        """
        if self._current_start is None:
            raise RuntimeError("Trying to close an already closed run of Ns.")
        if position <= self._current_start:
            raise RuntimeError(
                f"Trying to close a sequence of Ns with position {position} and start {self._current_start}"
            )
        self.runs.append((self._current_start, position - self._current_start))
        self._current_start = None


class Stats:
    def __init__(self, long_threshold: int, buckets_number, model: str) -> None:
        self._long_threshold = long_threshold
        self._buckets_number = buckets_number
        self._model = model

    def get_nbin(self, Sequence: dict):
        lines = [
            f"#Model {self._model} with >{self._long_threshold}bp of N Sequence\n",
            f"#SN\tBinID\tStart\tSize\n",
        ]
        for key, Sequence in Sequence.items():
            for index, run in enumerate(Sequence):
                lines.append(f"{key}\t{index+1}\t{run[0]:,}\t{run[1]:,}\n")
        return lines

    def get_nbuc(self):
        lines = [
            f"#WGS Extract generated Sequence of N summary file\n",
            f"#Model {self.model} with >{self._long_threshold}bp Sequence of N and {self._buckets_number} buckets per sequence\n",
            f"#Seq\tNumBP\tNumNs\tNumNreg\tNregSizeMean\tNregSizeStdDev\tSmlNreg\tBuckSize\t\n",
        ]


class UnknownBasesStats:

    def __init__(
        self, path: Path, long_run_threshold: int = 300, buckets_number: int = 1000
    ):
        self._long_run_threshold = long_run_threshold
        self._buckets_number = buckets_number
        self._path = path

    def _process_file(self, fp: TextIO) -> Dict[str, list[Tuple[int, int]]]:
        sequences = dict()
        pattern = re.compile(r"N+")

        position = None
        current_sequence = None

        for line in fp:
            line = line.rstrip("\n")

            if len(line) == 0 or line[0] == "#":
                continue

            if line[0] == "+":
                raise RuntimeError(
                    f"Expected a FASTA reference model, got a FASTQ sequencer data."
                )

            # Start of a new sequence
            if line[0] == ">":
                if current_sequence is not None and current_sequence.is_open():
                    current_sequence.close(position)
                sequence_name = line.split()[0][1:]
                logging.info(f"Processing sequence {sequence_name}")
                if sequence_name in sequences:
                    raise RuntimeError(f"Found a duplicated sequence: {sequence_name}")
                if "1" in sequences:
                    #    pass
                    return sequences
                current_sequence = Sequence(sequence_name)
                sequences[sequence_name] = current_sequence.runs
                position = 0
                continue

            result = list(pattern.finditer(line))

            if len(result) == 0:
                # Next sequence found but there's still an open run: close
                if current_sequence.is_open():
                    current_sequence.close(position)

            for match in result:
                if match.start() == 0:
                    # If open, keep open: is continuing from the previous line.
                    if not current_sequence.is_open():
                        current_sequence.open(position)
                else:
                    if current_sequence.is_open():
                        current_sequence.close(position)
                    current_sequence.open(match.start() + position)

                if match.end() < len(line):
                    current_sequence.close(match.end() + position)
            position += len(line)

        # File was terminated but there's still an open run: close
        if current_sequence.is_open():
            current_sequence.close(position)
        return sequences

    def _count_unknown_bases(self):
        logging.info(f"Starting to process file: {self._path.name}")
        with gzip.open(self._path, "rt") as f:
            return self._process_file(f)

    def _compute_statistics(self, unknown_basis: Dict[str, List[Tuple[int, int]]]):
        stats = Stats(self._long_run_threshold, self._buckets_number, self._path.stem)
        long_Sequence = {
            x: [y for y in unknown_basis[x] if y[1] >= self._long_run_threshold]
            for x in unknown_basis
        }
        lines = stats.get_nbin(long_Sequence)
        with open("checazzo.csv", "wt") as f:
            f.writelines(lines)

    def get_stats(self):
        unknown_bases = self._count_unknown_bases()
        self._compute_statistics(unknown_bases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("unknown_base_stat.py")
    parser.add_argument(
        "reference",
        help="Indicate the path of the bgzip compressed reference genome",
        type=Path,
    )
    parser.add_argument("--verbose", help="Set logging level [1,5]", type=int, default=2)
    args = parser.parse_args()
    logging.getLogger().setLevel(args.verbose)
    
    start = time.time()
    stats = UnknownBasesStats(Path("reference\\genomes\\hs37d5.fa.gz"))
    # countn = UnknownBasesStats(args.reference)
    statistics = stats.get_stats()
    end = time.time()
    logging.info(f"Time {end-start}s")
