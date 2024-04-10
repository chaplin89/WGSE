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
    
    def get_file_size(self, url:str):
        curl = pycurl.Curl()
        curl.setopt(pycurl.URL, url)
        curl.setopt(pycurl.NOBODY, 1)
        curl.perform()
        return curl.getinfo(pycurl.CONTENT_LENGTH_DOWNLOAD)
        

    def download(self, genome: Genome, callback: any = False) -> None:
        try:
            genome.initial_name.unlink(True)
            with genome.initial_name.open("wb") as f:
                curl = pycurl.Curl()
                curl.setopt(pycurl.URL, genome.url)
                curl.setopt(pycurl.WRITEDATA, f)
                curl.setopt(pycurl.FOLLOWLOCATION, True)
                if callback is not False or callback is not None:
                    curl.setopt(pycurl.NOPROGRESS, False)

                # TODO: find an easy way (if it exists) to load ca-bundle
                # on every platform. Or remove this TODO as we don't really
                # care about people MITM-ing our reference genome download.
                curl.setopt(pycurl.SSL_VERIFYPEER, 0)
                curl.setopt(pycurl.SSL_VERIFYHOST, 0)

                if callback is not None and callback is not False:
                    curl.setopt(pycurl.PROGRESSFUNCTION, lambda tot,n,tot_up,n_up: callback(tot, n))
                curl.perform()
                curl.close()
                # Signal the download is done
                callback(None, None)
        except:
            if genome.initial_name.exists():
                genome.initial_name.unlink()
            raise
