#!/usr/bin/env python

"""
Count the repeat length by measuring distance. Repeat for R1 and R2. When both R1 and R2 are read through reads, R2 results overwrite R1

Dev notes:
1. Nextflow won't resume if utils changes.
2. When determining the length, not simply use "left_start - right_end" because the adopted method is more robust for fuzzy match of the pattern.
3. If allow mismatch 2, length determined from R2 will be systematically 12 nt longer because of the fuzzy match: GGAGGCGGCGGCGGCGGCGG CGGCGGTGGCGG
    Solution:
        1. Use finditer(overlapped = True) to find all of the matches
        2. Use the one with no fuzzy_counts
        3. If 2 does not exist, use the right-most of the left match and left-most of the right match
        4. This would shrink to the minimal length
        Ref: https://stackoverflow.com/questions/5616822/how-to-use-regex-to-find-all-overlapping-matches
4. When matching, used "s<=", unlike "e<=" in the classify_readthrough.py, which is intentional. Therefore, it is possible to observe "NA,NA" results.

"""

from Bio import SeqIO
import sys
import regex
import gzip
from mimetypes import guess_type
from functools import partial
import numpy as np
import csv
from utils import reverse_complement

r1 = sys.argv[1]
r2 = sys.argv[2]
outfile = sys.argv[3]
r1_flanking = sys.argv[4]
r2_flanking = sys.argv[5]
mismatch = sys.argv[6]

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r1_flanking_rc = reverse_complement(r1_flanking)
r2_flanking_rc = reverse_complement(r2_flanking)

dict_repeat_length_per_read  = {}

# search flanking seqs in both R1 and R2, if both flankings exist, calculate the repeat length
with _open(r1) as f1, _open(r2) as f2:
    m = str(mismatch)
    for record_r1, record_r2 in zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')):
        search_seq_r1 = str(record_r1.seq)
        search_seq_r2 = str(record_r2.seq)

        # r1_left  = regex.search("".join(["(", r1_flanking, "){s<=", m, "}"]), search_seq_r1)
        # r1_right = regex.search("".join(["(", r2_flanking_rc, "){s<=", m, "}"]), search_seq_r1)
        #
        # r2_left  = regex.search("".join(["(", r2_flanking, "){s<=", m, "}"]), search_seq_r2)
        # r2_right = regex.search("".join(["(", r1_flanking_rc, "){s<=", m, "}"]), search_seq_r2)

        r1_left = None
        tem = regex.finditer("".join(["(", r1_flanking, "){s<=", m, "}"]), search_seq_r1, overlapped = True)
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                r1_left = match
                break
            r1_left = match

        r1_right = None
        tem = list(regex.finditer("".join(["(", r2_flanking, "){s<=", m, "}"]), search_seq_r1, overlapped = True))[::-1]
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                r1_right = match
                break
            r1_right = match
        # r1_right = next(regex.finditer("".join(["(", r2_flanking_rc, "){s<=", m, "}"]), search_seq_r1, overlapped = True), None)

        r2_left  = None
        tem  = regex.finditer("".join(["(", r2_flanking_rc, "){s<=", m, "}"]), search_seq_r2, overlapped = True)
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                r2_left = match
                break
            r2_left = match
        r2_right = None
        tem = list(regex.finditer("".join(["(", r1_flanking_rc, "){s<=", m, "}"]), search_seq_r2, overlapped = True))[::-1]
        for match in tem:
            if match.fuzzy_counts == (0, 0, 0):
                r2_right = match
                break
            r2_right = match
        # r2_right = next(regex.finditer("".join(["(", r1_flanking_rc, "){s<=", m, "}"]), search_seq_r2, overlapped = True), None)

        record_name = record_r1.name
        assert record_name == record_r2.name, "R1 and R2 sequence name differ!"

        if r1_left and r1_right:
            left_match_start = r1_left.start()
            left_match_length = len(r1_left.captures()[0])
            right_match_end = r1_right.end()
            right_match_length = len(r1_right.captures()[0])
            repeat_length = right_match_end - left_match_start - left_match_length - right_match_length
            dict_repeat_length_per_read[record_name] = [repeat_length]
        else:
            dict_repeat_length_per_read[record_name] = [np.nan]

        if r2_left and r2_right:
            left_match_start = r2_left.start()
            left_match_length = len(r2_left.captures()[0])
            right_match_end = r2_right.end()
            right_match_length = len(r2_right.captures()[0])
            repeat_length = right_match_end - left_match_start - left_match_length - right_match_length
            dict_repeat_length_per_read[record_name].append(repeat_length)
        else:
            dict_repeat_length_per_read[record_name].append(np.nan)

# ouput repeat_length_per_read_default:
with open(outfile, "w", newline = "") as f:
    writer = csv.writer(f)
    for key, value in dict_repeat_length_per_read.items():
        writer.writerow([key] + value)
