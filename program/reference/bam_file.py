import enum
import logging
import pathlib
import subprocess
import typing
import sys
sys.path.insert(0, ".\\")

from .external import External
from .genome import Genome
from .fasta_dictionary import ChromosomeNameType, FastaDictionary, Gender, MTDNANameType, Sorting
from .genome_repository import GenomeRepository

logger = logging.getLogger(__name__)


class ReadType(enum.Enum):
    SingleEnd = 0
    PariedEnd = 1
    Unknown = 2


class MTDNAType(enum.Enum):
    rCRS = 0
    Yoruba = 1
    RSRS = 2


class AlignmentMapFileType(enum.Enum):
    BAM = 0
    SAM = 1
    CRAM = 2


class ChromosomeType(enum.Enum):
    Autosomal = 0
    X = 1
    Y = 2
    Mitocondrial = 3
    Other = 4
    Unmapped = 5


class AlignmentMapFileInfo:
    def __init__(self) -> None:
        self.sorted: Sorting = None
        self.indexed: bool = None
        self.file_type: AlignmentMapFileType = None
        self.reference_genome: Genome = None
        self.read_type: ReadType = None
        self.is_y_only: bool = None
        self.is_mt_only: bool = None
        self.reference_mtdna: MTDNAType = None
        self.build: int = None
        self.name_type_chromosomes: ChromosomeNameType = None
        self.name_type_mtdna: MTDNANameType = None
        self.sequence_count: int = None
        self.primary: bool = None
        self.sequencer = None
        self.gender: Gender = None


class SequenceNameNormalizer:
    def normalize(name: str):
        return name.upper().replace("CHR", "").replace("MT", "M")


class SequenceStatistics:
    def __init__(self, type, name, reference_length, n_count, mapped) -> None:
        self.type = type
        self.name = name
        self.reference_length = reference_length
        self.n_count = n_count
        self.mapped = mapped


class AlignmentDepthStats:
    def __init__(self) -> None:
        pass


class AlignmentSummaryStats:
    def __init__(self, file: pathlib.Path, external: External = None) -> None:
        if not file.exists():
            raise RuntimeError(f"Unable to find file {file.name}")
        if external is None:
            external = External()

        self._file = file
        self._external = external
        #self.stats_by_type, self.lenght_by_type = self._get_stats()

    def _get_stats(self):
        stats_by_type: typing.Dict[ChromosomeType, SequenceStatistics] = dict()
        length_by_type: typing.Dict[ChromosomeType, int] = dict()

        stats_text = self._external.samtools(
            ["idxstats", self._file], stdout=subprocess.PIPE, wait=True
        )
        stats_text: list[str] = stats_text.decode().split("\n")
        for line in stats_text:
            # Line format is: name, length, # mapped, # unmapped
            chromosome_type = None
            line = line.strip()
            if line == "":
                continue
            elements = line.split("\t")
            name = elements[0]
            reference_length = int(elements[1])
            mapped = int(elements[2])
            unmapped = int(elements[3])
            read_lenght = mapped + unmapped

            name_normalized = SequenceNameNormalizer.normalize(name)

            if name_normalized.isnumeric():
                chromosome_type = ChromosomeType.Autosomal
            elif name_normalized == "M":
                chromosome_type = ChromosomeType.Mitocondrial
            elif name_normalized == "X":
                chromosome_type = ChromosomeType.X
            elif name_normalized == "Y":
                chromosome_type = ChromosomeType.Y
            elif name_normalized == "*":
                chromosome_type = ChromosomeType.Unmapped
            else:
                chromosome_type = ChromosomeType.Other

            stats = SequenceStatistics(
                chromosome_type, name, reference_length, 0, mapped
            )
            if chromosome_type not in stats_by_type:
                stats_by_type[chromosome_type] = list()
            stats_by_type[chromosome_type].append(stats)

            if chromosome_type not in length_by_type:
                length_by_type[chromosome_type] = 0
            length_by_type[chromosome_type] += read_lenght
            return stats_by_type, length_by_type

    def get_gender(self):
        x_length = self.lenght_by_type.get(ChromosomeType.X, 0)
        y_length = self.lenght_by_type.get(ChromosomeType.Y, 0)

        if x_length == 0 and y_length == 0:
            return None
        elif y_length == 0 or (x_length / y_length) > 20:
            return Gender.Female
        else:
            return Gender.Male


class AlignmentMapFile:
    SUPPORTED_FILES = {
        ".bam": AlignmentMapFileType.BAM,
        ".sam": AlignmentMapFileType.SAM,
        ".cram": AlignmentMapFileType.CRAM,
    }

    def __init__(self, file: pathlib.Path, external: External = None) -> None:
        if not file.exists():
            raise FileNotFoundError(f"Unable to find file {file.name}")
        if file.suffix.lower() not in AlignmentMapFile.SUPPORTED_FILES:
            raise RuntimeError(f"Unrecognized file extension: {file.name}")
        if external is None:
            external = External()

        self._file = file
        self._external = external
        self._header = self._load_header()
        self._summary = AlignmentSummaryStats(file, self._external)

        self._file_info = self._initialize_file_info()
        
    def _load_header(self) -> FastaDictionary:
        lines = self._external.samtools(
            ["view", "-H", "--no-PG", self._file], stdout=subprocess.PIPE, wait=True
        )
        lines = lines.decode().split("\n")
        dictionary = FastaDictionary(lines)
        return dictionary
        

    def _initialize_file_info(self):
        file_info = AlignmentMapFileInfo()
        file_info.file_type = AlignmentMapFile.SUPPORTED_FILES[
            self._file.suffix.lower()
        ]
        file_info.sorted = self._header.header.sorted
        file_info.name_type_mtdna = self._header.mtdna_name_type()
        file_info.name_type_chromosomes = self._header.chromosome_name_type()
        file_info.sequence_count = self._header.sequence_count()
        file_info.indexed = self.is_indexed(file_info.file_type)
        # file_info.gender = self._summary.get_gender()
        return file_info

    def is_indexed(self, type=None):
        if type == None:
            type = self._file_info.file_type

        file = str(self._file)
        if type == AlignmentMapFileType.BAM:
            return pathlib.Path(file + ".bai").exists()
        elif type == AlignmentMapFileType.CRAM:
            return pathlib.Path(file + ".crai").exists()
        return None


class DeterministicSequenceOrderer:
    def __init__(self, sequences=typing.List[str]) -> None:
        self._sequences = {self._normalize(x): x for x in sequences}
        self._deterministic_ordered = self._get_ordered()
        
    def __iter__(self):
        for sequence in self._deterministic_ordered:
            yield self._sequences[sequence]

    def _get_ordered(self):
        autosomal = self._get_autosomal()
        sexual = self._get_sexual()
        mitochondrial = self._get_mitocondrial()
        others = self._get_others(autosomal, sexual, mitochondrial)
        merged = [*autosomal, *sexual, *mitochondrial, *others]
        return merged

    def _get_autosomal(self):
        return [x for x in self._sequences if x.isnumeric()]
    
    def _get_mitocondrial(self):
        return [x for x in self._sequences if x == "m"]

    def _get_sexual(self):
        sexual = []
        if "x" in self._sequences:
            sexual.append("x")
        if "y" in self._sequences:
            sexual.append("y")
        return sexual

    def _get_others(self, autosomal, sexual, mitochondrial ):
        others = []
        for sequence in self._sequences.keys():
            is_autosomal = sequence in autosomal
            is_sexual = sequence in sexual
            is_mitochondial = sequence in mitochondrial
            if not is_autosomal and not is_sexual and not is_mitochondial:
                others.append(sequence)
        others.sort()
        return others

    def _normalize(self, sequence_name: str):
        normalized = sequence_name.lower()
        if normalized.startswith("chr"):
            normalized = normalized.replace("chr", "", 1)
        if normalized.startswith("mt"):
            normalized = normalized.replace("mt", "m", 1)
        return normalized


class ReferenceFinder:
    def __init__(self, repository: GenomeRepository = None) -> None:
        if repository is None:
            repository = GenomeRepository.build()
        self._repository = repository

    def find(self, map: AlignmentMapFile):
        genomes = self._repository.filter()
        valid_genomes = []
        for genome in genomes:
            if not genome.dict.exists():
                continue
            dictionary : FastaDictionary = FastaDictionary.load_from_file(genome.dict)
            valid_genomes.append(dictionary)
            sequences = DeterministicSequenceOrderer(list(dictionary.entries.keys()))
            am_sequences = DeterministicSequenceOrderer(map._header.entries.keys())
            
            ref_lenghts = [dictionary.entries[x].length for x in sequences]
            am_lenght = [map._header.entries[x].length for x in am_sequences]
            
            if ref_lenghts == am_lenght:
                print(f"Matches!!! {genome.code}")
            else:
                print(f"Not matching!!!11 {genome.code}")