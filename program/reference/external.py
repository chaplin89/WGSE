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

    def get_bio_default_directory():
        if "win" in sys.platform:
            return pathlib.Path("cygwin64", "usr", "local", "bin")
        else:
            return pathlib.Path("/", "usr", "bin")
    
    def get_sys_default_directory():
        if "win" in sys.platform:
            return pathlib.Path("cygwin64", "bin")
        else:
            return pathlib.Path("/", "bin")
    
    def __init__(self, installation_directory: Path = None) -> None:
        if installation_directory == None:
            installation_directory_bio = External.get_bio_default_directory()
            installation_directory_sys = External.get_sys_default_directory()
        else:
            installation_directory_sys = installation_directory_bio = installation_directory
        
        if not installation_directory_bio.exists():
            raise FileNotFoundError(
                f"Unable to find root directory for External: {str(installation_directory_bio)}"
            )
            
        self._installation_directory_bio = installation_directory_bio
        self._installation_directory_sys = installation_directory_sys
        self._htsfile = str(self._installation_directory_bio.joinpath("htsfile"))
        self._samtools = str(self._installation_directory_bio.joinpath("samtools"))
        self._tabix = str(self._installation_directory_bio.joinpath("tabix"))
        self._bgzip = str(self._installation_directory_bio.joinpath("bgzip"))
        self._bcftools = str(self._installation_directory_bio.joinpath("bcftools"))
        self._gzip = str(self._installation_directory_sys.joinpath("gzip"))

        files_collection = [
            self._htsfile,
            self._samtools,
            self._tabix,
            self._bgzip,
            self._bcftools,
            self._gzip
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
    
    def _gzip_filename(self, input:Path, action: GzipAction):
        if action == GzipAction.Compress:
            return Path(str(input) + ".gz")
        elif action == GzipAction.Decompress:
            if len(input.suffixes) == 0:
                raise RuntimeError(
                    f"Unable to determine decompressed filename, invalid filename {str(input)} (no extensions)."
                )
            return Path(str(input).rstrip(input.suffixes[-1]))
        elif action == GzipAction.Reindex:
            return Path(str(input) + ".gzi")
        else:
            raise RuntimeError(f"Action {action.name} not supported.")
            
    def gzip(
        self, input: Path, output: Path, action: GzipAction = GzipAction.Decompress) -> Path:
        if output.exists():
            raise RuntimeError(f"Trying to decompress {str(input)} but the destination file {str(output)} exists.")
        inferred_filename = self._gzip_filename(input, action)
        
        action_flags = {
            GzipAction.Compress:"",
            GzipAction.Decompress:"-d"
        }
        
        arguments = [self._gzip, action_flags[action], str(input)]
        process = subprocess.run(arguments, capture_output=True)
        
        if process.returncode != 0:
            if "trailing garbage" not in process.stderr.decode():
                raise RuntimeError(f"gzip failed: {process.stderr}")
        
        if inferred_filename != output:
            inferred_filename.rename(output)

    def bgzip(
        self, input: Path, output: Path, action: GzipAction = GzipAction.Compress, index: Path = None
    ) -> Path:
        if output.exists():
            raise RuntimeError(f"Trying to decompress {str(input)} but the destination file {str(output)} exists.")
        
        action_flags = {
            GzipAction.Compress:"-if",
            GzipAction.Decompress:"-d",
            GzipAction.Reindex: "-r"
        }
        inferred_filename = self._gzip_filename(input, action)

        arguments = [self._bgzip, action_flags[action], str(input), "-@", "32"]
        process = subprocess.run(arguments, capture_output=True)
        if inferred_filename != output:
            inferred_filename.rename(output)

    def tab_indexer(
        self,
    ):
        """Generic indexer for TAB-delimited genome position files"""
        arguments = [self._tabix]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")

    def idxstats(self, input: Path):
        """Generate BAM index statistics"""
        arguments = [self._samtools, "idxstat", input]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")
