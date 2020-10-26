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

def count_repeats(seq):
    """
    Convert sequence into pair of (seq, counts) where every nucleotide
    is paired with a count of how many times it repeats.
    """
    chars = []
    counts = []
    last_char = None
    count = 0
    for c in seq:
        if c == last_char:
            count += 1
        else:
            if last_char and count > 0:
                chars.append(last_char)
                counts.append(count)
            last_char = c
            count = 1
    return "".join(chars), np.array(counts)

