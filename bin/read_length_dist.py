#!/usr/bin/env python

"""
To output read length statistics:
file1:
    read_name, read_length

file2:
    read_length, number_of_reads

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
import numpy as np
import pandas as pd
from utils import plot_repeat_dist

fastq = sys.argv[1]
sample_name = sys.argv[2]

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

dict_count  = {}

with _open(fastq) as f:
    for record in SeqIO.parse(f, 'fastq'):
        dict_count[record.name] = str(record.seq)

for x in dict_count:
    f.write(x + "," + str(dict_count[x]) + "\n")
