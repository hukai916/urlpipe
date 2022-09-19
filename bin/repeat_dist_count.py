#!/usr/bin/env python

"""
To count the repeat length by measuring distance.

Usage:
    python count_cag.py r1.fastq dir_output

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
from collections import Counter
import regex
import gzip
from mimetypes import guess_type
from functools import partial
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import sys
import pandas as pd

r1 = sys.argv[1]
r2 = sys.argv[2]
sample_name = sys.argv[3]
output_dir  = sys.argv[4]
r1_flanking = sys.argv[5]
r2_flanking = sys.argv[6]
mismatch = sys.argv[7]

# r1_flanking = "TCGAGTCCCTCAAGTCCTTC"
# r2_flanking = "CCGCCACCGCCGCCGCCGCC"
# mismatch = 2
encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r1_dna = Seq(r1_flanking)
r2_dna = Seq(r2_flanking)
r1_flanking_rc = str(r1_dna.reverse_complement())
r2_flanking_rc = str(r2_dna.reverse_complement())

dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        r1_search = regex.search("(" + r1_flanking + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        r2_search = regex.search("(" + r2_flanking_rc + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        if not r1_search:
            dict_repeat["problem"] += 1 # should not happen
            dict_count[record.name] = "problem"
        elif not r2_search:
            dict_repeat["plus"] += 1
            dict_count[record.name] = "plus"
        else:
            r1_match_start = r1_search.start()
            r1_match_length = len(r1_search.captures()[0])
            r2_match_end = r2_search.end()
            r2_match_length = len(r2_search.captures()[0])
            repeat_length = r2_match_end - r1_match_start - r1_match_length - r2_match_length
            dict_count[record.name] = repeat_length

            if not repeat_length in dict_repeat:
                dict_repeat[repeat_length] = 1
            else:
                dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r1", sample_name + ".tsv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "\t" + str(dict_repeat[i]) + "\n")
    f.write("plus" + "\t" + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "\t" + str(dict_repeat["problem"]) + "\n")

# output plot:
output_plot = os.path.join(output_dir, "plot_r1", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# We can set the number of bins with the *bins* keyword argument.
df = pd.read_csv(output_file, sep = "\t")
x, y = [], []
_dict = {}

for i in range(len(df.iloc[:,0])):
    _dict[df.iloc[i, 0]] = df.iloc[i, 1]

for i in list(range(int(df.iloc[-3, 0]) + 10)) + ["plus", "problem"]:
    if not str(i) in _dict:
        x = x + [str(i)]
        y = y + [0]
    else:
        x = x + [str(i)]
        y = y + [_dict[str(i)]]

# weight = [math.log2(x + 0.1) for x in df.iloc[:, 1]]
weight = y
# n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
n, bins, patches = plt.hist(x, weights = weight, bins=len(x), edgecolor='black')

# sparse the x-axis ticks:
N = 5  # 1 tick every 5
myticks = [i for i in range(len(x)) if not i%N]
newlabels = [i for i in range(len(x)) if not i%N]

if "plus" in list(x) and not "plus" in newlabels:
    newlabels = newlabels[:-2]
    myticks = myticks[:-2]
    newlabels = newlabels + ["+"]
    myticks = myticks + [len(x) - 2]

plt.xticks(myticks, newlabels, rotation = 85)
plt.tick_params(labelsize=8)
plt.title(sample_name)

plt.savefig(output_plot, dpi = 300)
plt.close()

# ouput raw count:
output_count = os.path.join(output_dir, "count_r1", sample_name + ".tsv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "\t" + str(dict_count[x]) + "\n")


## Now use R2:
dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        r2_search = regex.search("(" + r2_flanking + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        r1_search = regex.search("(" + r1_flanking_rc + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        if not r2_search:
            dict_repeat["problem"] += 1 # should not happen
            dict_count[record.name] = "problem"
        elif not r1_search:
            dict_repeat["plus"] += 1
            dict_count[record.name] = "plus"
        else:
            r2_match_start = r2_search.start()
            r2_match_length = len(r2_search.captures()[0])
            r1_match_end = r1_search.end()
            r1_match_length = len(r1_search.captures()[0])
            repeat_length = r1_match_end - r2_match_start - r2_match_length - r1_match_length
            dict_count[record.name] = repeat_length

            if not repeat_length in dict_repeat:
                dict_repeat[repeat_length] = 1
            else:
                dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r2", sample_name + ".tsv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "\t" + str(dict_repeat[i]) + "\n")
    f.write("plus" + "\t" + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "\t" + str(dict_repeat["problem"]) + "\n")

# output plot:
output_plot = os.path.join(output_dir, "plot_r2", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# We can set the number of bins with the *bins* keyword argument.
df = pd.read_csv(output_file, sep = "\t")
x, y = [], []
_dict = {}

for i in range(len(df.iloc[:,0])):
    _dict[df.iloc[i, 0]] = df.iloc[i, 1]

for i in list(range(int(df.iloc[-3, 0]) + 10)) + ["plus", "problem"]:
    if not str(i) in _dict:
        x = x + [str(i)]
        y = y + [0]
    else:
        x = x + [str(i)]
        y = y + [_dict[str(i)]]

# weight = [math.log2(x + 0.1) for x in df.iloc[:, 1]]
weight = y
# n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
n, bins, patches = plt.hist(x, weights = weight, bins=len(x), edgecolor='black')

# sparse the x-axis ticks:
N = 5  # 1 tick every 5
myticks = [i for i in range(len(x)) if not i%N]
newlabels = [i for i in range(len(x)) if not i%N]

if "plus" in list(x) and not "plus" in newlabels:
    newlabels = newlabels[:-2]
    myticks = myticks[:-2]
    newlabels = newlabels + ["+"]
    myticks = myticks + [len(x) - 2]

plt.xticks(myticks, newlabels, rotation = 85)
plt.tick_params(labelsize=8)
plt.title(sample_name)

plt.savefig(output_plot, dpi = 300)
plt.close()

# ouput raw count:
output_count = os.path.join(output_dir, "count_r2", sample_name + ".tsv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "\t" + str(dict_count[x]) + "\n")
