from .fasta_dictionary import FastaDictionary
import collections
import gzip
import logging
import re
import typing

from .genome import Genome

class Run:
    """Represent a single run of Ns"""

    def __init__(self) -> None:
        self.start = None
        self.length = None

    def open(self, position: int):
        self.start = position

    def close(self, length: int):
        self.length = length

class Sequence:
    """Represent a collection of runs of Ns."""

    def __init__(self, name: str, length: int) -> None:
        if length <= 0:
            raise IndexError("Length should be greater than zero.")
        self.runs: typing.List[Run] = []
        self.name: str = name
        self.length: int = length
        self._current_run: Run = None
        self._is_ended: bool = False

    def is_run_open(self) -> bool:
        return self._current_run is not None

    def open_run(self, position: int) -> None:
        if self._is_ended:
            raise RuntimeError("Trying to open a run on an ended sequence.")
        if self._current_run is not None:
            raise RuntimeError("Trying to open a an already opened run of Ns.")
        if position < 0:
            raise IndexError("Trying to open a run with a negative position.")

        if len(self.runs) > 0:
            previous_end = self.runs[-1].start + self.runs[-1].length
            if position < previous_end + 1:
                raise ValueError("Trying to open a run overlapping with previous runs.")

        self._current_run = Run()
        self._current_run.open(position)

    def filter(self, criteria: typing.Callable[[Run], bool]):
        return [x for x in self.runs if criteria(x)]

    def close_run(self, position: int) -> None:
        if self._current_run is None:
            raise RuntimeError("Trying to close an already closed run of Ns.")
        if position <= self._current_run.start:
            raise RuntimeError(
                f"Expected a position greater than start for closing a run: {position}<={self._current_run.start}"
            )
        if position > self.length:
            raise RuntimeError(
                f"Trying to close a run with a position {position}, that is greater than the whole length of the sequence {self.length}."
            )
        run_length = position - self._current_run.start
        self._current_run.close(run_length)
        self.runs.append(self._current_run)
        self._current_run = None

    def end(self, position: int) -> None:
        if self.is_run_open():
            self.close_run(position)
        self._is_ended = True

        if position != self.length:
            raise ValueError(
                f"Expected {self.length} base pairs in this sequence but processed {position}."
            )

class FastaFile:
    """Represent a collection of sequences"""

    def __init__(
        self,
        genome: Genome
    ):
        self.genome = genome
        if not self.genome.dict.exists():
            raise RuntimeError(
                f"Unable to find dictionary in {self.genome.dict.name}."
            )
        self._dict = FastaDictionary(self.genome.dict)

    @property
    def model_name(self):
        return self.genome.final_name.name

    def _sequence_from_line(self, line: str) -> Sequence:
        # Processing the opening of a new sequence in a FASTA file.
        # It's usually a line that begins with '>' followed by the name
        # of the sequence.
        sequence_name = line.split()[0][1:]
        if sequence_name not in self._dict.entries:
            raise ValueError(f"Sequence {sequence_name} is not present in dictionary.")
        sequence_length = self._dict.entries[sequence_name].length
        return Sequence(sequence_name, sequence_length)

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

            # Check if processing a .fastq file by mistake.
            if line[0] == "+":
                raise RuntimeError(
                    "Expected a FASTA reference model, got a FASTQ sequencer data."
                )

            if line[0] == ">":
                # New sequence found. Close old sequence if open.
                if current_sequence is not None:
                    current_sequence.end(position)

                current_sequence = self._sequence_from_line(line)
                if current_sequence.name in sequences:
                    raise RuntimeError(
                        f"Found a duplicated sequence: {current_sequence.name}"
                    )
                logging.debug(f"{self.genome.final_name.name}: Processing sequence {current_sequence.name}")
                sequences[current_sequence.name] = current_sequence
                position = 0
                continue

            matches = list(pattern.finditer(line))

            # No Ns found in this line: close the open run (if any).
            if len(matches) == 0:
                if current_sequence.is_run_open():
                    current_sequence.close_run(position)

            for match in matches:
                if match.start() == 0:
                    # The Ns starts at the beginning of the line.
                    # If a run is open, keep open: is continuing from the previous line.
                    # Otherwise open a new one.
                    if not current_sequence.is_run_open():
                        current_sequence.open_run(position)
                else:
                    # There are Ns but they don't start at the beginning of the line.
                    # Close the run (if open) and re-open it: this is a new run of Ns.
                    if current_sequence.is_run_open():
                        current_sequence.close_run(position)
                    current_sequence.open_run(match.start() + position)

                # The run is ending before the end of the line. Close the open run (if any).
                if match.end() < len(line):
                    current_sequence.close_run(match.end() + position)
            position += len(line)

        # File is terminated: close the current sequence and the open run (if any).
        if current_sequence is None:
            raise RuntimeError("Found the end of the file but no sequences found.")
        current_sequence.end(position)
        return list(sequences.values())

    def split_into_sequences(self) -> typing.List[Sequence]:
        logging.info(f"{self.genome.final_name.name}: Counting Ns.")
        with gzip.open(self.genome.final_name, "rt") as f:
            sequences = self._process_file(f)
        logging.info(f"{self.genome.final_name.name}: Counting Ns done.")
        return sequences