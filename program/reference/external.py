import enum
from pathlib import Path
import pathlib
import subprocess
import sys


class GzipAction(enum.Enum):
    Compress = 0
    Decompress = 1
    Reindex = 2


class External:
    """Wrapper around External files"""

    def get_default_directory():
        if "win" in sys.platform:
            return pathlib.Path("cygwin64", "usr", "local", "bin")
        else:
            return pathlib.Path("usr", "bin")
    
    def __init__(self, installation_directory: Path = None) -> None:
        if installation_directory == None:
            installation_directory = External.get_default_directory()
        
        if not installation_directory.exists():
            raise FileNotFoundError(
                f"Unable to find root directory for External: {str(installation_directory)}"
            )
            
        self._installation_directory = installation_directory
        self._htsfile = str(self._installation_directory.joinpath("htsfile"))
        self._samtools = str(self._installation_directory.joinpath("samtools"))
        self._tabix = str(self._installation_directory.joinpath("tabix"))
        self._bgzip = str(self._installation_directory.joinpath("bgzip"))
        self._bcftools = str(self._installation_directory.joinpath("bcftools"))
        self._gunzip = str(self._installation_directory.joinpath("gzip"))

        files_collection = [
            self._htsfile,
            self._samtools,
            self._tabix,
            self._bgzip,
            self._bcftools,
        ]
        unix_files = all([x for x in files_collection if Path(x).exists()])
        win_files = all([x for x in files_collection if Path(x + ".exe").exists()])
        if not (unix_files or win_files):
            raise FileNotFoundError("Unable to find all the required 3rd party tools.")

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
            output = Path(path.parent, path.name + ".dict")
        arguments = [self._samtools, "dict", str(path), "-o", str(output)]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def gzip(
        self, path: Path, action: GzipAction = GzipAction.Compress, index: Path = None
    ) -> Path:
        """Wrapper around bgzip executable, used to compress a file
        in the bgzip format.

        Args:
            path (Path): Path for the file to be compressed.
            output (Path, optional): Output file path. Defaults to None.

        Returns:
            str: Output of the bgzip command.
        """
        action_flag = None
        if action == GzipAction.Compress:
            action_flag = ""
            output = Path(str(path) + ".gz")
        elif action == GzipAction.Decompress:
            action_flag = "-d"
            if len(path.suffixes) == 0:
                raise RuntimeError(
                    f"Unable to decompress, invalid filename {str(path)}"
                )

            output = Path(path.name.rstrip(path.suffixes[-1]))
        else:
            raise RuntimeError(f"Action {action.name} not supported.")

        arguments = [self._gunzip, action_flag, str(path)]
        process = subprocess.run(arguments, capture_output=True)
        return output

    def bgzip(
        self, path: Path, action: GzipAction = GzipAction.Compress, index: Path = None
    ) -> Path:
        """Wrapper around bgzip executable, used to compress a file
        in the bgzip format.

        Args:
            path (Path): Path for the file to be compressed.
            output (Path, optional): Output file path. Defaults to None.

        Returns:
            str: Output of the bgzip command.
        """
        action_flag = None
        if action == GzipAction.Compress:
            action_flag = "-if"
            output = [Path(str(path) + ".gz"), Path(str(path) + ".gzi")]
        elif action == GzipAction.Decompress:
            action_flag = "-d"
            if len(path.suffixes) == 0:
                raise RuntimeError(
                    f"Unable to decompress, invalid filename {str(path)}"
                )

            output = Path(path.name.rstrip(path.suffixes[-1]))
        elif action == GzipAction.Reindex:
            action_flag = "-r"
            output = Path(str(path) + ".gzi")

        arguments = [self._bgzip, action_flag, str(path), "-@", "32"]
        process = subprocess.run(arguments, capture_output=True)
        return output

    def tab_indexer(
        self,
    ):
        """Generic indexer for TAB-delimited genome position files"""
        arguments = [self._tabix]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")

    def idxstats(self, input: Path):
        """Generate BAM index statistics"""
        arguments = [self._External, input]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")
