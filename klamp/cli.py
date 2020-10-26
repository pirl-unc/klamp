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


import sys

import argparse
from glob import glob

from .dna import reverse_complement
from .fastq import FastQ
from .fasta import read_virus_reference
from .kmer import KmerIndex
from .repeats import count_repeats
from .version import __version__

parser = argparse.ArgumentParser("klamp")

parser.add_argument(
    "--sample",
    nargs="*",
    help="FASTQ file(s) containing ONT reads from LAMP amplification",
    required=True)

parser.add_argument(
    "--reference",
    help="FASTA file containing viral genome",
    required=True)

def expand_sample_paths(sample_paths):
    if type(sample_paths) is str:
        sample_paths = [sample_paths]
    result = []
    for sample_path in sample_paths:
        if "*" in sample_path:
            result.extend(glob(sample_path))
        else:
            result.append(sample_path)
    return result

def run_from_args(args, verbose=False):
    reference = read_virus_reference(args.reference)
    print("Read %0.1fkb reference from %s" % (
        len(reference) / 1000, args.reference,))

    compressed_reference, compressed_reference_counts = count_repeats(reference)
    print("Collapsed repeats in reference: %0.1fkb" % (
        (len(compressed_reference)/1000),
    ))
    compressed_revcomp = reverse_complement(compressed_reference)

    kmer_index = KmerIndex()
    kmer_index.index(compressed_reference)
    kmer_index.index(compressed_revcomp)
    sample_paths = expand_sample_paths(args.sample)
    sample_dict = {}
    for sample_path in sample_paths:
        print("Reading sample '%s'..." % sample_path,)
        fastq = FastQ.read(sample_path)
        sample_dict[sample_path] = fastq
        print("Read %d lines (%0.1fkb) from sample '%s'" % (
            fastq.num_lines(),
            fastq.num_bases() / 1000,
            sample_path))
        sample_dict[sample_path] = fastq

    for sample_path in sample_paths:
        fastq = sample_dict[sample_path]
        n_hits_total = 0
        n_kmers_total = 0
        n_reads = 0
        for i, line in enumerate(fastq.sequences()):
            n_reads += 1
            compressed_line, _ = count_repeats(line)

            if len(compressed_line) < kmer_index.kmer_size:
                continue

            n_hits, n_kmers = kmer_index.count_hits(compressed_line)
            n_hits_total += n_hits
            n_kmers_total += n_kmers

            if n_hits > 0:
                if verbose:
                    print("Line #%d: homopolymer compression %d => %d bases" % (
                        i + 1,
                        len(line),
                        len(compressed_line),
                    ))
                    print("%s:%10d - %d/%d hits (%0.2f%%)" % (
                        sample_path,
                        i + 1,
                        n_hits,
                        n_kmers,
                        n_hits / n_kmers * 100,
                    ))
        print("%s (reads=%d), %d/%d kmers match reference: %0.2f%%" % (
            sample_path,
            n_reads,
            n_hits_total,
            n_kmers_total,
            n_hits_total / n_kmers_total * 100))


def main(args_list=None):
    print("klamp version %s" % __version__)
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)
    run_from_args(args)

