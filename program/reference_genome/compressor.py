from genome import Genome
from samtools import Samtools
from pathlib import Path


class Compressor:

    def __init__(self, samtools: Samtools) -> None:
        self._samtools = samtools

    def need_compression(self, genome: Genome):
        pass
    
    def _is_bgzip_compressed(self, file: Path):
        """Uses samtools to check if a file is bgzip compressed.
        
        Args:
            file (Path): File to check

        Returns:
            bool: True if bgzip compressed, False otherwise.
        """
        file_type = self._samtools.get_file_type(file)

        if "FASTA BGZF" in file_type:
            return True
        return False

    def compress(self, genome: Genome):
        try:
            self._samtools.bgzip_compress(genome.uncompressed, genome.bgzip_compressed)
            
        except:
            if genome.bgzip_compressed.exists():
                genome.bgzip_compressed.unlink()