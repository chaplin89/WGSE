from .fasta_file import Sequence

import collections
import math
import typing


class Buckets:
    """Represent a sequence partitioned into subsequences of fixed lenght called buckets"""

    def __init__(
        self, sequence: Sequence, buckets_number: int, long_run_threshold: int
    ) -> None:
        self._buckets_number = buckets_number
        self._long_run_threshold: int = long_run_threshold
        self._sequence: Sequence = sequence
        self.buckets = self._make_buckets()

    def _make_buckets(self) -> typing.OrderedDict[int, int]:
        if self._buckets_number > self._sequence.length:
            # Can't have less than 1 N per bucket.
            return collections.OrderedDict()

        buckets = collections.OrderedDict()
        bucket_size = int(math.floor(self._sequence.length / self._buckets_number))

        for run in self._sequence.filter(
            lambda x: x.length >= self._long_run_threshold
        ):
            # Determine how many buckets this run is spanning
            start = run.start
            end = run.start + run.length
            index_bucket_start = int(math.floor(start / bucket_size))
            index_bucket_end = int(math.floor(end / bucket_size))

            # Iterate over each bucket determining how many
            for bucket in range(index_bucket_start, index_bucket_end + 1):
                start_bucket_offset = max(bucket * bucket_size, start)
                end_bucket_offset = min(bucket * bucket_size + bucket_size - 1, end)
                runs_count = end_bucket_offset - start_bucket_offset
                if runs_count < 1:
                    continue
                if bucket not in buckets:
                    buckets[bucket] = 0
                buckets[bucket] += runs_count
        return buckets

    def __getitem__(self, key):
        return self.buckets[key]