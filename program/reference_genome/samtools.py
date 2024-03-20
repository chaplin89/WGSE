from pathlib import Path
import subprocess


class Samtools:
    """Wrapper around Samtools files"""

    def __init__(self, installation_directory: Path) -> None:
        self._installation_directory = installation_directory
        self._htsfile = self._installation_directory.joinpath("htsfile")
        self._samtools = self._installation_directory.joinpath("samtools")
        self._tabix = self._installation_directory.joinpath("tabix")
        self._bgzip = self._installation_directory.joinpath("bgzip")
        self._bcftools = self._installation_directory.joinpath("bcftools")

    def get_file_type(self, path: Path):
        process = subprocess.run([self._htsfile, path], capture_output=True, check=True)
        return process.stdout.decode("utf-8")

    def fasta_index(self, path: Path, output: Path = None):
        """Create an index for a FASTA file.

        Args:
            path (Path): Path of the FASTA file.
            output (Path, optional): Target output file. Defaults to None.
        """
        if output is None:
            output = Path(str(path) + ".fai")

        arguments = [self._samtools, "faidx", path, "-o", output]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def view(self, file: Path, output: Path):
        arguments = [self._samtools, "view", "-H", "--no-PG", file]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def make_dictionary(self, path: Path, output: Path = None):
        """Wrapper around "samtool dict" command. Create a sequence dictionary
        file from a fasta file.

        Args:
            fasta_file (Path): Path of the FASTA (eventually bgzip compressed) file.
            output (Path, optional): Target output file. Defaults to None.

        Returns:
            str: Standard output of the samtool command.
        """
        if output is None:
            output = Path(str(path) + ".dict")
        arguments = [self._samtools, "dict", path, "-o", output]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def bgzip_compress(self, path: Path, output: Path = None):
        """Wrapper around bgzip executable, used to compress a file
        in the bgzip format.

        Args:
            path (Path): Path for the file to be compressed.
            output (Path, optional): Output file path. Defaults to None.

        Returns:
            str: Output of the bgzip command.
        """
        arguments = [self._bgzip, "-if", path]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def tab_indexer(self,):
        """Generic indexer for TAB-delimited genome position files"""
        arguments = [self._tabix]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")

    def idxstats(self, input:Path):
        """Generate BAM index statistics"""
        arguments = [self._samtools, input]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")