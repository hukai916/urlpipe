#!/usr/bin/env python

"""
To annotate reads on_target, off_target, and problem
    on_target: both ends match
    off_target: neither ends match
    problem (mis_primed): one end matches

Usage:
    python map_locus.py r1.fastq r2.fastq on_target_dir off_target_dir problem_dir stat_dir r1_flanking r2_flanking
"""

from Bio import SeqIO
import sys
import os
from collections import Counter
import gzip
from mimetypes import guess_type
from functools import partial
import regex

r1 = sys.argv[1]
r2 = sys.argv[2]
on_target_dir = sys.argv[3]
off_target_dir = sys.argv[4]
problem_dir = sys.argv[5]
stat_dir    = sys.argv[6]
r1_flanking = sys.argv[7]
r2_flanking = sys.argv[8]
error       = sys.argv[9]

r1_match = {}
r2_match = {}
r_match  = {}

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r1_flanking + ")" + "{e<=" + str(error) + "}", str(record.seq[:len(r1_flanking)])):
            r1_match[record.name] = 1
        else:
            r1_match[record.name] = 0
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r2_flanking + ")" + "{e<=" + str(error) + "}", str(record.seq[:len(r2_flanking)])):
            r2_match[record.name] = 1
        else:
            r2_match[record.name] = 0

for k in r1_match:
    assert k in r2_match, "not equal records!"
    r_match[k] = r1_match[k] + r2_match[k]

# output:
out_on_target_r1  = os.path.join(on_target_dir, os.path.basename(r1).split(".gz")[0])
out_on_target_r2  = os.path.join(on_targe_dir, os.path.basename(r2).split(".gz")[0])
out_off_target_r1 = os.path.join(off_target_dir, os.path.basename(r1).split(".gz")[0])
out_off_target_r2 = os.path.join(off_target_dir, os.path.basename(r2).split(".gz")[0])
out_problem_r1 = os.path.join(problem_dir, os.path.basename(r1).split(".gz")[0])
out_problem_r2 = os.path.join(problem_dir, os.path.basename(r2).split(".gz")[0])

for dir in [out_on_target_r1, out_on_target_r2, out_off_target_r1, out_off_target_r2, out_problem_r1, out_problem_r2]:
    os.makedirs(os.path.dirname(dir), exist_ok=True)

r1_on_target = []
r1_off_target = []
r1_problem = []
with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r1_on_target.append(record)
        elif r_match[record.name] == 1:
            r1_problem.append(record)
        elif r_match[record.name] == 0:
            r1_off_target.append(record)

SeqIO.write(r1_on_target, out_on_target_r1, "fastq")
SeqIO.write(r1_off_target, out_on_target_r1, "fastq")
SeqIO.write(r1_problem, out_problem_r1, "fastq")

r2_on_target = []
r2_off_target = []
r2_problem = []
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r2_on_target.append(record)
        elif r_match[record.name] == 1:
            r2_problem.append(record)
        elif r_match[record.name] == 0:
            r2_off_target.append(record)
SeqIO.write(r2_on_target, out_on_target_r2, "fastq")
SeqIO.write(r2_off_target, out_on_target_r2, "fastq")
SeqIO.write(r2_problem, out_problem_r2, "fastq")

# print some stats:
res = Counter(r_match.values())

os.makedirs(os.path.dirname(stat_dir), exist_ok = True)
with open(stat_dir, "w") as f:
    name = os.path.basename(r1).split("_1.fastq.gz")[0]
    p2 = str(res[2]/(sum([res[2], res[1], res[0]])))
    p1 = str(res[1]/(sum([res[2], res[1], res[0]])))
    p0 = str(res[0]/(sum([res[2], res[1], res[0]])))

    f.write(name + "," + str(res[2]) + "," + p2 + "," + str(res[1]) + "," + p1 + "," + str(res[0]) + "," + p0 + "\n")
