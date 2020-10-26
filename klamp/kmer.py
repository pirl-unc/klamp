# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import Counter

class KmerIndex(object):
    def __init__(self, kmer_size=25):
        self.kmer_size = kmer_size
        self.kmer_counts = Counter()

    def index(self, seq):
        n = len(seq)
        k = self.kmer_size
        counts = self.kmer_counts
        for i in range(n - k + 1):
            counts[seq[i:i + k]] += 1

    def count_hits(self, seq, stride=None):
        """
        Returns pair of integers:
        - number of kmers matching index
        - total number of kmers extracted
        """
        if stride is None:
            stride = self.kmer_size // 2
        n = len(seq)
        k = self.kmer_size
        index_kmers = self.kmer_counts
        n_kmers = 0
        n_hits = 0
        for i in range(0, n - k + 1, stride):
            kmer = seq[i:i + k]
            assert len(kmer) == k
            n_kmers += 1
            if kmer in index_kmers:
                n_hits += 1

        return n_hits, n_kmers

    def frac_hits(self, seq):
        """
        Returns fraction of kmers which match the index in this sequence
        """
        n_hits, n_total = self.count_hits(seq)
        return n_hits / n_total