#!/usr/bin/env python

"""
To plot the UMI distribution pattern.

"""

import sys
import os
from collections import Counter
from statistics import mean
from utils import plot_repeat_dist
from Bio import SeqIO
from Bio.Seq import Seq

import regex
import gzip
from mimetypes import guess_type
from functools import partial
import numpy as np
import pandas as pd
from utils import plot_repeat_dist
import copy


csv         = sys.argv[1]
sample_name = sys.argv[2]
read1       = sys.argv[3]
read2       = sys.argv[4]
outdir      = sys.argv[5]
cutoffs     = sys.argv[6]

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

umi_dict = {}
# per cutoff:
with open(csv, "r") as f:
    for line in f:
        umi, count, mean, mode, ld = line.split(",")
        count, mean, mode, ld = int(count), float(mean), int(mode), float(ld)
        umi_dict[umi] = [count, mean, mode, ld]

def process_cutoff(sample_name, cutoff, umi_dict):
    outfile_fastq_r1_mode = os.path.join(outdir, cutoff, "mode", sample_name, "_1.fastq.gz")
    outfile_fastq_r2_mode = os.path.join(outdir, cutoff, "mode", sample_name, "_2.fastq.gz")
    outfile_count_mode = os.path.join(outdir, cutoff, "mode", sample_name, ".csv")

    outfile_fastq_r1_ld = os.path.join(outdir, cutoff, "ld", sample_name, "_1.fastq.gz")
    outfile_fastq_r2_ld = os.path.join(outdir, cutoff, "ld", sample_name, "_2.fastq.gz")
    outfile_count_ld = os.path.join(outdir, cutoff, "ld", sample_name, ".csv")

    os.makedirs(os.path.join(outdir, cutoff, "mode"), exist_ok = True)
    os.makedirs(os.path.join(outdir, cutoff, "ld"), exist_ok = True)

    record_ld_r1 = []
    record_ld_r2 = []
    umi_dict_work = copy.deepcopy(umi_dict)
    with _open(read1) as f:
        for record in SeqIO.parse(f, 'fastq'):
            umi = record.name.split("_")[-1]
            if umi in umi_dict_work:
                if umi_dict_work[umi][0] >= cutoff and umi_dict_work[umi][3] == len(str(record.seq)):
                    record_ld_r1.append(record)
                    del umi_dict_work[umi]
    umi_dict_work = copy.deepcopy(umi_dict)
    with _open(read2) as f:
        for record in SeqIO.parse(f, 'fastq'):
            umi = record.name.split("_")[-1]
            if umi in umi_dict_work:
                if umi_dict_work[umi][0] >= cutoff and umi_dict_work[umi][3] == len(str(record.seq)):
                    record_ld_r2.append(record)
                    del umi_dict_work[umi]
    SeqIO.write(record_ld_r1, outfile_fastq_r1_ld, "fastq")
    SeqIO.write(record_ld_r2, outfile_fastq_r2_ld, "fastq")
    with open(outfile_count_ld, "w") as f:
        f.write(sample_name, ",", len(outfile_fastq_r1_ld))

    record_mode_r1 = []
    record_mode_r2 = []
    umi_dict_work = copy.deepcopy(umi_dict)
    with _open(read1) as f:
        for record in SeqIO.parse(f, 'fastq'):
            umi = record.name.split("_")[-1]
            if umi in umi_dict_work:
                if umi_dict_work[umi][0] >= cutoff and umi_dict_work[umi][2] == len(str(record.seq)):
                    record_ld_r1.append(record)
                    del umi_dict_work[umi]
    umi_dict_work = copy.deepcopy(umi_dict)
    with _open(read2) as f:
        for record in SeqIO.parse(f, 'fastq'):
            umi = record.name.split("_")[-1]
            if umi in umi_dict_work:
                if umi_dict_work[umi][0] >= cutoff and umi_dict_work[umi][2] == len(str(record.seq)):
                    record_ld_r2.append(record)
                    del umi_dict_work[umi]
    SeqIO.write(record_mode_r1, outfile_fastq_r1_mode, "fastq")
    SeqIO.write(record_mode_r2, outfile_fastq_r2_mode, "fastq")
    with open(outfile_count_mode, "w") as f:
        f.write(sample_name, ",", len(outfile_fastq_r1_mode))

for cutoff in cutoffs.split(","):
    process_cutoff(sample_name, int(cutoff), umi_dict)
