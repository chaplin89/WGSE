import logging
import zipfile
import gzip
import shutil
from pathlib import Path
from .file_type_checker import Type, FileTypeChecker


class Decompressor:
    def __init__(self, type_checker: FileTypeChecker) -> None:
        self._type_checker = type_checker

        self._handlers = {
            Type.GZIP: Decompressor.gz,
            Type.ZIP: Decompressor.zip,
            Type.SEVENZIP: Decompressor.sevenzip,
            Type.BZIP: Decompressor.bzip,
        }

    def gz(self, input_file: Path, output_file: Path):
        logging.debug(f"Decompressing file {input_file.name}. gzip compression detected.")

        with gzip.open(str(input_file), "rb") as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def sevenzip(self, input_file: Path, final_name: Path):
        logging.debug(f"Decompressing file {input_file.name}. 7z compression detected.")
        raise NotImplementedError()

    def bzip(self, input_file: Path, output_file: Path):
        logging.debug(f"Decompressing file {input_file}. bzip compression detected.")
        raise NotImplementedError()

    def zip(self, input_file: Path, output_file: Path):
        logging.debug(f"Decompressing file {input_file.name}. zip compression detected.")
        with zipfile.ZipFile(str(input_file), "r") as f:
            files = f.namelist()
            if len(files) > 1:
                raise RuntimeError("More than one file found inside the .zip, unable to proceed.")
            extracted = Path(f.extract(files[0], output_file.parent))
            if extracted != output_file:
                extracted.rename(output_file)

    def _get_target_file(self, input_file: Path):
        output_file = Path(str(input_file) + "_decompressed")
        return output_file

    def decompress(self, input_file: Path):
        """Decompress a genome file.

        Args:
            genome (Genome): Genome to decompress.

        Raises:
            FileNotFoundError: Input compressed genome file is not found.
            RuntimeError: Extension is not recognized.
        """
        output_file = self._get_target_file(input_file)
        if output_file.exists():
            output_file.unlink()
        
        if not input_file.exists():
            raise FileNotFoundError(
                f"Unable to find reference genome file {input_file.name} to decompress."
            )
        type = self._type_checker.get_type(input_file)
        if type not in self._handlers:
            raise RuntimeError(
                f"Trying to decompress a reference genome file with an unknown file extension: {input_file.name}"
            )

        handler = self._handlers[type]
        handler(self, input_file, output_file)
        if not output_file.exists():
            raise RuntimeError(f"Unable to decompress file {input_file.name} into {output_file.name}.")
        return output_file