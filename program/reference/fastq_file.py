import os
import re
from math import sqrt

from utilities import DEBUG, nativeOS
import settings as wgse


class FASTQFile:
    """
    Central handler for FASTQ files. Processes, determines stats, and manipulates data.

    Attributes:
        filename (str): Path to the FASTQ file.
        paired (bool): Whether the FASTQ file is paired-end.
        seqid (str): Sequencer ID identified from the file.
        num_segments (int): Number of segments in the file.
        avg_read_length (float): Average read length in the file.
        avg_read_stddev (float): Standard deviation of read lengths.
    """

    def __init__(self, filename, paired=True):
        """
        Initializes a FASTQFile object.

        Args:
            filename (str): Path to the FASTQ file.
            paired (bool, optional): Whether the FASTQ file is paired-end. Defaults to True.
        """

        self.filename = nativeOS(filename)
        self.paired = paired
        self.seqid = None
        self.num_segments = None
        self.avg_read_length = None
        self.avg_read_stddev = None

        self._process_fastq()

    def _process_fastq(self):
        """
        Processes the FASTQ file and calculates statistics.
        """

        line_cnt = char_cnt = 0
        rcnt = rlen = rmean = rM2 = 0

        with gzip.open(self.filename, 'rt') as f:
            for line in f:
                if line_cnt % 4 == 2:
                    rlen += len(line)
                    if rlen > 1:
                        rcnt += 1
                        rdelta = rlen - rmean
                        rmean += rdelta / rcnt
                        rdelta2 = rlen - rmean
                        rM2 += rdelta * rdelta2
                elif line_cnt == 0:
                    if line[0] == '#':
                        continue
                    else:
                        self.seqid = determine_sequencer(re.split("[ \t]", line.strip())[0][1:])
                        line_cnt = 1
                elif line_cnt > 20000:
                    break
                char_cnt += len(line)
                line_cnt += 1

        if rcnt > 2:
            rstd = sqrt(rM2 / (rcnt - 1))
            self.avg_read_length = rmean
            self.avg_read_stddev = rstd

        fastq_stats = os.stat(self.filename)
        if fastq_stats.st_size > 0:
            self.num_segments = int(fastq_stats.st_size / (char_cnt / (line_cnt / 4)))
            self.num_segments *= 2 if self.paired else 1

        DEBUG(
            f'FASTQ Stats: ID - "{self.seqid}, # segs - {self.num_segments:,d},'
            f' avg read length - {self.avg_read_length:,.0f}, read len stddev - {self.avg_read_stddev:,.0f}')


def determine_sequencer(seqid):
    """
    Determines sequencer type from the FASTQ SeqID.

    Args:
        seqid (str): The sequence ID.

    Returns:
        str: The sequencer type or "Error" if not found.
    """

    if seqid is None:
        return "Error"

    for key, val in wgse.sequencers.items():
        if re.compile(val[0]).match(seqid):
            if key != "Unknown":
                return key
    return seqid[0:23]+"..." if len(seqid) > 25 else seqid[0:25]


# Example usage
fastq_file = FASTQFile("path/to/your.fastq")
print(f"Sequencer ID: {fastq_file.seqid}")
print(f"Number of segments: {fastq_file.num_segments}")
print(f"Average read length: {fastq_file.avg_read_length}")
