import logging
import typing
import zipfile
import gzip
import shutil
from pathlib import Path
from .external import External, GzipAction
from .file_type_checker import Type, FileTypeChecker


class Decompressor:
    def __init__(self, type_checker: FileTypeChecker, external: External) -> None:
        self._type_checker = type_checker

        self._handlers : typing.Dict[Type, typing.Callable[[Path, Path], None]] = {
            Type.GZIP: Decompressor.gz,
            Type.ZIP: Decompressor.zip,
            Type.SEVENZIP: Decompressor.sevenzip,
            Type.BZIP: Decompressor.bzip,
            Type.RAZF_GZIP: Decompressor.razf_gzip
        }
        
        self._external = external

    def gz(self, input_file: Path, output_file: Path):
        logging.debug(f"Decompressing file {input_file.name}. gzip compression detected.")

        with gzip.open(str(input_file), "rb") as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def sevenzip(self, input_file: Path, output_file: Path):
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
    
    def razf_gzip(self, input_file: Path, output_file: Path):
        logging.debug(f"Decompressing file {input_file.name}. razf (gzip) compression detected.")
        self._external.gzip(input_file, output_file, GzipAction.Decompress)

    def decompress(self, input_file: Path, output_file: Path):
        if output_file.exists():
            output_file.unlink()
        
        if not input_file.exists():
            raise FileNotFoundError(
                f"Unable to find file {input_file.name} to decompress."
            )
        type = self._type_checker.get_type(input_file)
        if type not in self._handlers:
            raise RuntimeError(
                f"Trying to decompress a file with an unknown file extension: {input_file.name}"
            )

        handler = self._handlers[type]
        handler(self, input_file, output_file)
        if not output_file.exists():
            raise RuntimeError(f"Unable to decompress file {input_file.name} into {output_file.name}.")