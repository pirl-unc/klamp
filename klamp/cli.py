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

from .version import __version__
from .fastq import FastQ
from .fasta import read_virus_reference
from .repeats import count_repeats

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

def run_from_args(args):
    reference = read_virus_reference(args.reference)
    print("Read %0.1fkb reference from %s" % (
        len(reference) / 1000, args.reference,))

    compressed_reference, compressed_reference_counts = count_repeats(reference)
    print("Collapsed repeats in reference: %0.1fkb" % (
        (len(compressed_reference)/1000),
    ))
    sample_paths = args.sample
    if type(sample_paths) is str:
        sample_paths = [sample_paths]
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


def main(args_list=None):
    print("klamp version %s" % __version__)
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)
    run_from_args(args)

