import enum
from pathlib import Path
import pathlib
import subprocess
import sys


class GzipAction(enum.Enum):
    Compress = 0
    Decompress = 1
    Reindex = 2


def bin(f):
    def execute_binary(self, args=[], decode=True):
        if not isinstance(args, list):
            # Handle (common) case of a single parameter.
            args = [str(args)]
        
        process = External.get_bin_folder().joinpath(f.__name__)
        args = [process, *[str(x) for x in args]]
        output = subprocess.run(args, capture_output=True)
        if output.returncode != 0:
            raise RuntimeError(
                f"Process {f.__name__} failed with error: {output.stderr}"
            )
        if decode:
            return output.stdout.decode()
        return output.stdout
    return execute_binary

def usr_bin(f):
    def execute_binary(self, args=[], decode=True):
        if not isinstance(args, list):
            # Handle (common) case of a single parameter.
            args = [args]
        
        process = External.get_usr_bin_folder().joinpath(f.__name__)
        args = [process, *args]
        output = subprocess.run(args, capture_output=True)
        if output.returncode != 0:
            raise RuntimeError(
                f"Process {f.__name__} failed with error: {output.stderr}"
            )
        if decode:
            return output.stdout.decode()
        return output.stdout
    return execute_binary

def bio(f):
    def execute_binary(self, args=[], decode=True):
        if not isinstance(args, list):
            # Handle (common) case of a single parameter.
            args = [args]
        
        process = External.get_bio_folder().joinpath(f.__name__)
        args = [process, *args]
        output = subprocess.run(args, capture_output=True)
        if output.returncode != 0:
            raise RuntimeError(
                f"Process {f.__name__} failed with error: {output.stderr}"
            )
        if decode:
            return output.stdout.decode()
        return output.stdout
    return execute_binary

class External:
    """Wrapper around External files"""

    def get_bio_folder():
        if "darwin" in sys.platform.lower():
            return pathlib.Path("/", "opt", "local", "bin")
        return External.get_usr_bin_folder()

    def get_usr_bin_folder():
        if "win" in sys.platform:
            return pathlib.Path("cygwin64", "usr", "local", "bin")
        return pathlib.Path("/", "usr", "bin")

    def get_bin_folder():
        if "win" in sys.platform:
            return pathlib.Path("cygwin64", "bin")
        return pathlib.Path("/", "bin")

    def get_java8_folder():
        # Windows
        jre8 = f'{install_FP}jre8/bin/java'     # If Win10 installer had to install, then local to WGSE
        if not java8x_FN and is_command_available(jre8, "-version", True):
            java8x_FN = jre8
            java8_version = cp_version
        jre17 = f'{install_FP}jre17/bin/java'   # If Win10 installer had to install, then local to WGSE
        if not java17x_FN and is_command_available(jre17, "--version", True):
            java17x_FN = jre17
            java17_version = cp_version
        # Mac
        jre8 = f'/Library/Java/JavaVirtualMachines/zulu-8.jre/Contents/Home/bin/java'
        if not java8x_FN and is_command_available(jre8, "-version", True):
            java8x_FN = jre8
            java8_version = cp_version
        
        if not java17x_FN:
            import glob
            jvms = glob.glob(f'/Library/Java/JavaVirtualMachines/zulu-[12]?.jre/Contents/Home')
            if len(jvms) > 0:
                jre17 = f'{jvms[0]}/bin/java'
                if is_command_available(jre17, "--version", True):
                    java17x_FN = jre17
                    java17_version = cp_version
        # Linux
        jre8 = f'/usr/lib/jvm/java-8-openjdk-amd64/bin/java'
        if not java8x_FN and is_command_available(jre8, "-version", True):
            java8x_FN = jre8
            java8_version = cp_version
        if not java17x_FN:
            import glob
            jvms = glob.glob(f'/usr/lib/jvm/java-[12]?-openjdk-amd64')
            if len(jvms) > 0:
                jre17 = f'{jvms[0]}/bin/java'
                if is_command_available(jre17, "--version", True):
                    java17x_FN = jre17
                    java17_version = cp_version

    def __init__(self, installation_directory: Path = None) -> None:
        if installation_directory == None:
            usr_local = External.get_usr_bin_folder()
            bin = External.get_bin_folder()
        else:
            bin = usr_local = installation_directory

        if not usr_local.exists():
            raise FileNotFoundError(
                f"Unable to find root directory for External: {str(usr_local)}"
            )

        self._usr_bin = usr_local
        self._bin = bin
        self._htsfile = str(self._usr_bin.joinpath("htsfile"))
        self._samtools = str(self._usr_bin.joinpath("samtools"))
        self._bgzip = str(self._usr_bin.joinpath("bgzip"))
        self._gzip = str(self._bin.joinpath("gzip"))
        
        files_collection = [
            self._htsfile,
            self._samtools,
            self._bgzip,
        ]

        unix_files = all([x for x in files_collection if Path(x).exists()])
        win_files = all([x for x in files_collection if Path(x + ".exe").exists()])
        if not (unix_files or win_files):
            raise FileNotFoundError("Unable to find all the required 3rd party tools.")

    def get_file_type(self, path: Path):
        process = subprocess.run([self._htsfile, path], capture_output=True, check=True)
        return process.stdout.decode("utf-8")

    def fasta_index(self, path: Path, output: Path = None):
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

    @bio
    def bwa(self,*k):
        pass
    @bio
    def bwamem2(self,*k):
        pass
    @bio
    def minimap2(self,*k):
        pass
    @bio
    def fastp(self,*k):
        pass
    @bio
    def bcftool(self, *k):
        pass
    @bio
    def tabix(self, *k):
        pass

    @usr_bin
    def head(self, *k):
        pass
    @usr_bin
    def tail(self, *k):
        pass
    @usr_bin
    def gawk(self, *k):
        pass
    @usr_bin
    def grep(self, *k):
        pass
    @usr_bin
    def sort(self, *k):
        pass
    
    @bin
    def cat(self, *k):
        pass
    @bin
    def zcat(self, *k):
        pass
    @bin
    def wc(self, *k):
        pass
    @bin
    def sed(self, *k):
        pass
    @bin
    def zip(self, *k):
        pass
    @bin
    def unzip(self, *k):
        pass
    @bin
    def mv(self, *k):
        pass
    @bin
    def cut(self, *k):
        pass
    @bin
    def uniq(self, *k):
        pass
    @bin
    def pr(self, *k):
        pass