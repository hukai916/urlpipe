#!/usr/bin/env python


from Bio import SeqIO
import sys
import os
import gzip
from mimetypes import guess_type
from functools import partial
from utils import plot_repeat_dist

r1 = sys.argv[1]
r2 = sys.argv[2]
sample_name = sys.argv[3]
output_dir  = sys.argv[4]
bin_number = sys.argv[5]

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        repeat_length = len(str(record.seq))
        dict_count[record.name] = repeat_length

        if not repeat_length in dict_repeat:
            dict_repeat[repeat_length] = 1
        else:
            dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r1", sample_name + ".stat.csv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "," + str(dict_repeat[i]) + "\n")
    f.write("plus" + "," + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "," + str(dict_repeat["problem"]) + "\n")

# output plot:
N = 5

output_plot = os.path.join(output_dir, "plot_r1", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plot_repeat_dist(output_file, output_plot, sample_name, N, bin_number)

# ouput raw count:
output_count = os.path.join(output_dir, "count_r1", sample_name + ".csv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "," + str(dict_count[x]) + "\n")

## Now use R2:
dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        repeat_length = len(str(record.seq))
        dict_count[record.name] = repeat_length

        if not repeat_length in dict_repeat:
            dict_repeat[repeat_length] = 1
        else:
            dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r2", sample_name + ".stat.csv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "," + str(dict_repeat[i]) + "\n")
    f.write("plus" + "," + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "," + str(dict_repeat["problem"]) + "\n")

# output plot:
output_plot = os.path.join(output_dir, "plot_r2", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plot_repeat_dist(output_file, output_plot, sample_name, N, bin_number)

# ouput raw count:
output_count = os.path.join(output_dir, "count_r2", sample_name + ".csv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "," + str(dict_count[x]) + "\n")
