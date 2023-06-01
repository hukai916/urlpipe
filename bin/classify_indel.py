#!/usr/bin/env python

"""
To annotate reads into no_indel, indel_5p, indel_3p, and indel_5p_3p (flanking)
    # align R1 to 20bp before CAG repeat, allow 1 mismatch: TCGAGTCCCTCAAGTCCTTC
    # align R2 to 20bp after CAG repeat, allow 1 mismatch: CCGCCACCGCCGCCGCCGCC (use rev comp)
    # if both match: "no indel"; if R1 doesn't match: "indel, 5p"; if R2 doesn't match: "indel, 3p"
    # implementation: instead of perform mapping, simple search sub_string with python regex:
    regex.search("(xxgyy){s<=1}", "xxggyy") # https://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string, note, regex.search("xxgyy{s<=1}", "xxggyy") does not work.
Usage:
    classify_indel.py ${prefix}_1.fastq.gz ${prefix}_2s.fastq.gz no_indel indel_5p indel_3p indel_5p_3p stat sample_name $args

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
from utils import indel_filter

r1 = sys.argv[1]
r2 = sys.argv[2]
no_indel_dir = sys.argv[3]
indel_5p_dir = sys.argv[4]
indel_3p_dir = sys.argv[5]
indel_5p_and_3p_dir = sys.argv[6]
indel_5p_or_3p_dir = sys.argv[7]
indel_stat_dir  = sys.argv[8]
sample_name = sys.argv[9]
r1_flanking = sys.argv[10]
r2_flanking = sys.argv[11]
mismatch    = sys.argv[12]
indel_cutoff = float(sys.argv[13])

r1_match = {}
r2_match = {}
r_match  = {}

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r2_dna = Seq(r2_flanking)
r2_flanking_rc = str(r2_dna.reverse_complement())

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r1_flanking + ")" + "{s<=" + str(mismatch) + "}", str(record.seq)):
            r1_match[record.name] = 1
        else:
            r1_match[record.name] = 0
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r2_flanking_rc + ")" + "{s<=" + str(mismatch) + "}", str(record.seq)):
            r2_match[record.name] = 1
        else:
            r2_match[record.name] = 0

for k in r1_match:
    assert k in r2_match, "not equal records!"
    r_match[k] = r1_match[k] + r2_match[k]

# output:
output_dirs = {
    "no_indel": no_indel_dir,
    "indel_5p": indel_5p_dir,
    "indel_3p": indel_3p_dir,
    "indel_5p_and_3p": indel_5p_and_3p_dir,
    "indel_5p_or_3p": indel_5p_or_3p_dir
}
output_files = {}

for output_type, output_dir in output_dirs.items():
    output_files[output_type] = {
        "r1": os.path.join(output_dir, sample_name + "_1.fastq"),
        "r2": os.path.join(output_dir, sample_name + "_2.fastq")
    }

for outfile in output_files:
    os.makedirs(os.path.dirname(output_files[outfile]["r1"]), exist_ok = True)
    os.makedirs(os.path.dirname(output_files[outfile]["r2"]), exist_ok = True)

count_5p = 0
count_3p = 0

r1_no_indel = []
r1_indel_5p = []
r1_indel_3p = []
r1_indel_5p_and_3p = []
r1_indel_5p_or_3p = []
with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r1_no_indel.append(record)
        elif r_match[record.name] == 1:
            r1_indel_5p_or_3p.append(record)
            if r1_match[record.name]:
                r1_indel_3p.append(record)
                count_3p += 1
            else:
                r1_indel_5p.append(record)
                count_5p += 1
        elif r_match[record.name] == 0:
            r1_indel_5p_and_3p.append(record)
            r1_indel_5p_or_3p.append(record)

r1_indel_3p_filter = indel_filter(r1_indel_3p, r1_no_indel)
r1_indel_5p_filter = indel_filter(r1_indel_5p, r1_no_indel)
r1_indel_5p_and_3p_filter = indel_filter(r1_indel_5p_and_3p, r1_no_indel)
r1_indel_5p_or_3p_filter, r1_no_indel_filter = indel_filter(r1_indel_5p_or_3p, r1_no_indel, add = True)

_res = [r1_no_indel_filter, r1_indel_5p_filter, r1_indel_3p_filter, r1_indel_5p_and_3p_filter, r1_indel_5p_or_3p_filter]
_out_file = [output_files["no_indel"], output_files["indel_5p"], output_files["indel_3p"], output_files["indel_5p_and_3p"], output_files["indel_5p_or_3p"]]

for res, out_file in zip(_res, _out_file):
    with open(str(out_file["r1"]), "w") as f:
        SeqIO.write(res, f, "fastq")

r2_no_indel = []
r2_indel_5p = []
r2_indel_3p = []
r2_indel_5p_and_3p = []
r2_indel_5p_or_3p = []
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r2_no_indel.append(record)
        elif r_match[record.name] == 1:
            r2_indel_5p_or_3p.append(record)
            if r2_match[record.name]:
                r2_indel_5p.append(record)
            else:
                r2_indel_3p.append(record)
        elif r_match[record.name] == 0:
            r2_indel_5p_and_3p.append(record)
            r2_indel_5p_or_3p.append(record)

r2_indel_3p_filter = indel_filter(r2_indel_3p, r2_no_indel, indel_cutoff = indel_cutoff)
r2_indel_5p_filter = indel_filter(r2_indel_5p, r2_no_indel, indel_cutoff = indel_cutoff)
r2_indel_5p_and_3p_filter = indel_filter(r2_indel_5p_and_3p, r2_no_indel, indel_cutoff = indel_cutoff)
r2_indel_5p_or_3p_filter, r2_no_indel_filter = indel_filter(r2_indel_5p_or_3p, r2_no_indel, add = True, indel_cutoff = indel_cutoff)

_res = [r2_no_indel_filter, r2_indel_5p_filter, r2_indel_3p_filter, r2_indel_5p_and_3p_filter, r2_indel_5p_or_3p_filter]

for res, out_file in zip(_res, _out_file):
    with open(str(out_file["r2"]), "w") as f:
        SeqIO.write(res, f, "fastq")

# print some stats:
res = Counter(r_match.values())
with open(os.path.join(indel_stat_dir, sample_name + ".csv"), "w") as f:
    p2 = str(res[2]/(sum([res[2], res[1], res[0]])))
    p_5p = str(count_5p/(sum([res[2], res[1], res[0]])))
    p_3p = str(count_3p/(sum([res[2], res[1], res[0]])))
    p0 = str(res[0]/(sum([res[2], res[1], res[0]])))
    res = ",".join([sample_name, str(res[2]), str(count_5p), str(count_3p), str(res[0]), p2, p_5p, p_3p]) + "\n"

    f.write(res)
