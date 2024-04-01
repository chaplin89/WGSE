import enum
from pathlib import Path
import pathlib
import subprocess
import sys
import os

if "win" in sys.platform:
    usr_bin = str(Path(".", "cygwin64", "usr", "local", "bin"))
    bin = str(Path(".", "cygwin64", "bin"))
    mingw = str(Path(".", "usr", "local", "mingw64.bin"))
    if usr_bin not in os.environ["PATH"]:
        os.environ["PATH"] += ";" + usr_bin
    if bin not in os.environ["PATH"]:
        os.environ["PATH"] += ";" + bin
    if mingw not in os.environ["PATH"]:
        os.environ["PATH"] += ";" + mingw

class GzipAction(enum.Enum):
    Compress = 0
    Decompress = 1
    Reindex = 2


def run(f):
    def execute_binary(self, args=[], stdout=None, stdin=None, wait=False):
        if not isinstance(args, list):
            # Handle (common) case of a single parameter.
            args = [str(args)]
        args = [f.__name__, *[str(x) for x in args]]
        
        output = subprocess.Popen(args, stdout=stdout, stdin=stdin, stderr=subprocess.PIPE)
        if wait == True:
            out, err = output.communicate()
            if output.returncode!= 0:
                raise RuntimeError(f"Call to {f.__name__} failed: {err.decode()}")
            return out
        return output
    return execute_binary

class External:
    """Wrapper around External files"""

    def __init__(self, installation_directory: Path = None) -> None:
        if installation_directory is not None:
            if not installation_directory.exists():
                raise FileNotFoundError(
                    f"Unable to find root directory for External: {str(installation_directory)}"
                )
            if str(installation_directory) not in os.environ["PATH"]:
                os.environ["PATH"] += ";" + str(installation_directory)

        self._htsfile = "htsfile"
        self._samtools = "samtools"
        self._bgzip = "bgzip"
        self._gzip = "gzip"

        # files_collection = [
        #     self._htsfile,
        #     self._samtools,
        #     self._bgzip,
        # ]
        # unix_files = all([x for x in files_collection if Path(x).exists()])
        # win_files = all([x for x in files_collection if Path(x + ".exe").exists()])
        # if not (unix_files or win_files):
        #     raise FileNotFoundError("Unable to find all the required 3rd party tools.")

    def get_file_type(self, path: Path):
        process = subprocess.run([self._htsfile, path], capture_output=True, check=True)
        return process.stdout.decode("utf-8")

    def fasta_index(self, path: Path, output: Path = None):
        if output is None:
            output = Path(str(path) + ".fai")

        arguments = [self._samtools, "faidx", path, "-o", output]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def view(self, file: Path, output: Path, *args):
        arguments = [self._samtools, "view", "-H", "--no-PG", *args, file]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def make_dictionary(self, path: Path, output: Path = None):
        if output is None:
            output = Path(path.parent, path.name + ".dict")
        arguments = [self._samtools, "dict", str(path), "-o", str(output)]
        process = subprocess.run(arguments, check=True, capture_output=True)
        return process.stdout.decode("utf-8")

    def _gzip_filename(self, input: Path, action: GzipAction):
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
        self, input: Path, output: Path, action: GzipAction = GzipAction.Decompress
    ) -> Path:
        if output.exists():
            raise RuntimeError(
                f"Trying to decompress {str(input)} but the destination file {str(output)} exists."
            )
        inferred_filename = self._gzip_filename(input, action)

        action_flags = {GzipAction.Compress: "", GzipAction.Decompress: "-d"}

        arguments = [self._gzip, action_flags[action], str(input)]
        process = subprocess.run(arguments, capture_output=True)

        if process.returncode != 0:
            # RAFZ format is libz compatible but will make gzip returning
            # with a != 0 code, complaining about "trailing garbage data".
            # This is not a real error, as the file is decompressed anyway.
            # The issue is potentially fixable by truncating the file, but
            # there's no practical advantage in doing so. If we fall in this
            # situation, ignore the error.
            if "trailing garbage" not in process.stderr.decode():
                raise RuntimeError(f"gzip failed: {process.stderr}")

        if inferred_filename != output:
            inferred_filename.rename(output)

    def bgzip(
        self, input: Path, output: Path, action: GzipAction = GzipAction.Compress
    ) -> Path:
        if output.exists():
            raise RuntimeError(
                f"Trying to decompress {str(input)} but the destination file {str(output)} exists."
            )

        action_flags = {
            GzipAction.Compress: "-if",
            GzipAction.Decompress: "-d",
            GzipAction.Reindex: "-r",
        }
        inferred_filename = self._gzip_filename(input, action)

        arguments = [self._bgzip, action_flags[action], str(input), "-@", "32"]
        process = subprocess.run(arguments, capture_output=True)
        if inferred_filename != output:
            inferred_filename.rename(output)

    def idxstats(self, input: Path):
        """Generate BAM index statistics"""
        arguments = [self._samtools, "idxstat", input]
        process = subprocess.run(arguments)
        return process.stdout.decode("utf-8")

    @run
    def samtools(self, *k):
        pass

    @run
    def bwa(self, *k):
        pass

    @run
    def bwamem2(self, *k):
        pass

    @run
    def minimap2(self, *k):
        pass

    @run
    def fastp(self, *k):
        pass

    @run
    def bcftool(self, *k):
        pass

    @run
    def tabix(self, *k):
        pass

    @run
    def head(self, *k):
        pass

    @run
    def tail(self, *k):
        pass

    @run
    def gawk(self, *k):
        pass

    @run
    def grep(self, *k):
        pass

    @run
    def sort(self, *k):
        pass

    @run
    def cat(self, *k):
        pass

    @run
    def zcat(self, *k):
        pass

    @run
    def wc(self, *k):
        pass

    @run
    def sed(self, *k):
        pass

    @run
    def zip(self, *k):
        pass

    @run
    def unzip(self, *k):
        pass

    @run
    def mv(self, *k):
        pass

    @run
    def cut(self, *k):
        pass

    @run
    def uniq(self, *k):
        pass

    @run
    def pr(self, *k):
        pass
