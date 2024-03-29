import pycurl
from .genome import Genome
import hashlib

try:
    import tqdm
except:
    tqdm = None

class Downloader:
    def __init__(self) -> None:
        self._progressbar = None

    def callback(download_total, downloaded, upload_total, uploaded):
        pass
    
    def need_download(self, genome: Genome) -> bool:
        """Return True if a Genome need to be downloaded.

        Args:
            genome (Genome): Genome to check.

        Returns:
            bool: True if the genome needs to be downloaded, False otherwise.
        """
        if not genome.initial_name.exists():
            hashlib.md5(open('filename.exe','rb').read()).hexdigest()
            return True
        return False

    def download(self, genome: Genome, callback: any = False) -> None:
        """Download a genome, overwrite if it exist.

        Args:
            genome (Genome): Genome to download
            callback (any, optional): Function to call for progress.
                None disable progress reporting, False enable pycurl
                default progress reporting. Defaults to False.
        """

        try:
            genome.initial_name.unlink(True)
            with genome.initial_name.open("wb") as f:
                curl = pycurl.Curl()
                curl.setopt(pycurl.URL, genome.url)
                curl.setopt(pycurl.WRITEDATA, f)
                curl.setopt(pycurl.FOLLOWLOCATION, True)
                if callback is not False:
                    curl.setopt(pycurl.NOPROGRESS, False)

                # TODO: find an easy way (if it exists) to load ca-bundle
                # on every platform. Or remove this TODO as we don't really
                # care about people MITM-ing our reference genome download.
                curl.setopt(pycurl.SSL_VERIFYPEER, 0)
                curl.setopt(pycurl.SSL_VERIFYHOST, 0)

                if callback is not None and callback is not False:
                    curl.setopt(pycurl.DEBUGFUNCTION, callback)
                curl.perform()
                curl.close()
        except:
            if genome.initial_name.exists():
                genome.initial_name.unlink()
            raise
