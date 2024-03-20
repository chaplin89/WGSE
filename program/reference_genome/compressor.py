from genome import Genome
from samtools import Samtools
from pathlib import Path

class Compressor:
    def __init__(self, samtools: Samtools) -> None:
        self._samtools = samtools

    def determine_target_name(self, file: Path):
        if file.suffix == ".gz":
            return file
        return Path(str(file) + ".gz")

    def compress(self, genome: Genome):
        try:
            target = self.determine_target_name(genome.initial_name)
            self._samtools.bgzip_compress(genome.initial_name, target)
        except:
            if target.exists():
                target.unlink()
            raise
        return target