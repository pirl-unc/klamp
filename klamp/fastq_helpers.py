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

def read_id_and_metadata_line(line, convert_int=True):
    parts = line.split()
    name = parts[0]
    if name.startswith("@"):
        name = name[1:]
    else:
        raise ValueError("Malformed FASTQ file")
    metadata = {}
    for part in parts[1:]:
        if part.count("=") == 1:
            k, v = part.split("=")
            if convert_int and v.isdigit():
                v = int(v)
            metadata[k] = v
    return name, metadata

def read_quals(line):
    n = len(line)
    quals = np.zeros(n, dtype='int')
    for i, c in enumerate(line):
        quals[i] = ord(c) - 33
    return quals