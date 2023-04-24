#!/usr/bin/env python

"""
Count the repeat length by measuring distance using the merged reads.


"""

from Bio import SeqIO
import sys
import os
from collections import Counter
import regex
import gzip
from mimetypes import guess_type
from functools import partial
import numpy as np
import pandas as pd
import csv
from utils import reverse_complement, plot_repeat_dist

read = sys.argv[1]
outfile = sys.argv[2]
r1_flanking = sys.argv[3]
r2_flanking = sys.argv[4]
mismatch = sys.argv[5]

encoding = guess_type(read)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r1_flanking_rc = reverse_complement(r1_flanking)
r2_flanking_rc = reverse_complement(r2_flanking)

dict_repeat_length_per_read  = {}

# search flanking seqs in read, if both flankings exist, calculate the repeat length
with _open(read) as f:
    m = str(mismatch)
    for record in SeqIO.parse(f, 'fastq'):
        search_seq = str(record.seq)

        read_left = None
        tem = regex.finditer("".join(["(", r1_flanking, "){s<=", m, "}"]), search_seq, overlapped = True)
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                read_left = match
                break
            read_left = match

        read_right = None
        tem = list(regex.finditer("".join(["(", r2_flanking_rc, "){s<=", m, "}"]), search_seq, overlapped = True))[::-1]
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                read_right = match
                break
            read_right = match

        record_name = record.name

        # when both R1 and R2 are read through reads, R2 results overwrite R1
        if read_left and read_right:
            left_match_start = read_left.start()
            left_match_length = len(read_left.captures()[0])
            right_match_end = read_right.end()
            right_match_length = len(read_right.captures()[0])
            repeat_length = right_match_end - left_match_start - left_match_length - right_match_length
            dict_repeat_length_per_read[record_name] = [repeat_length]
        else:
            dict_repeat_length_per_read[record_name] = [np.nan]

# ouput repeat_length_per_read_merge.csv:
with open(outfile, "w", newline = "") as f:
    writer = csv.writer(f)
    for key, value in dict_repeat_length_per_read.items():
        writer.writerow([key] + value)