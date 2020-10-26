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

from .common import normalize_string

class Buffer():
    def __init__(self):
        self.seq_dict = {}
        self.reset()

    def reset(self):
        self.curr_name = None
        self.lines = []

    def finalize_current_record(self):
        self.seq_dict[self.curr_name] = "".join(self.lines)
        self.reset()

    def start_new_record(self, name):
        self.finalize_current_record()
        self.curr_name = name

    def add_line(self, line):
        self.lines.append(line)


def read_fasta(path):
    buffer = Buffer()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                buffer.start_new_record(line[1:].split()[0])
            else:
                buffer.add_line(normalize_string(line))
    buffer.finalize_current_record()
    return buffer.seq_dict

def read_virus_reference(path, seq_name="SARSCoV2"):
    seq_dict = read_fasta(path)
    if len(seq_dict) == 1:
        return list(seq_dict.values())[0]
    if seq_name in seq_dict:
        return seq_dict[seq_name]
    # search for similar key
    for k, v in seq_dict.items():
        if normalize_string(k) == normalize_string(seq_name):
            return v
    raise ValueError(
        "Unable to find sequence with name similar to '%s' in file '%s'" % (
            seq_name,
            path))