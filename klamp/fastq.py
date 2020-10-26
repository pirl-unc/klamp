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

import numpy as np

from .common import normalize_string
from .fastq_helpers import read_id_and_metadata_line, read_quals

class FastQ(object):
    def __init__(self, seq_dict, qual_dict, metadata_dict):
        self.seq_dict = seq_dict
        self.qual_dict = qual_dict
        self.metadata_dict = metadata_dict

    def qualities(self):
        return self.qual_dict.values()

    def keys(self):
        return self.seq_dict.keys()

    def sequences(self):
        return self.seq_dict.values()

    def key_set(self):
        return set(self.keys())

    def num_lines(self):
        return len(self.seq_dict)

    def num_bases(self):
        total = 0
        for s in self.sequences():
            total += len(s)
        return total

    @classmethod
    def read(cls, path, skip_quals=True):
        seq_dict = {}
        qual_dict = {}
        metadata_dict = {}
        with open(path) as f:

            curr_name = None
            curr_metadata = None
            curr_seq = None
            curr_quals = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if curr_name is None:
                    curr_name, curr_metadata = \
                        read_id_and_metadata_line(line)
                elif curr_seq is None:
                    curr_seq = normalize_string(line)
                elif line == "+":
                    continue
                elif curr_quals is None:
                    if skip_quals:
                        curr_quals = np.zeros(len(curr_seq), dtype="int")
                    else:
                        curr_quals = read_quals(line)
                    if len(curr_quals) != len(curr_seq):
                        raise ValueError(
                            "Malformed FASTQ file, read and quality lengths don't match")
                    seq_dict[curr_name] = curr_seq
                    qual_dict[curr_name] = curr_quals
                    metadata_dict[curr_name] = curr_metadata
                    curr_name = curr_seq = curr_quals = curr_metadata = None
                else:
                    raise ValueError("Malformed FASTQ")
        return cls(
            seq_dict=seq_dict,
            qual_dict=qual_dict,
            metadata_dict=metadata_dict)

