#!/usr/bin/env python

"""
Output UMI corrected (at different cutoff) reads in fastq.gz format.

Usage:
python get_reads_umi_cutoff.py umi_cutoff read_count_${prefix}_umi_cutoff_\$i.csv ${prefix}_1.fastq.gz reads_per_umi_${prefix}.csv out1.fastq out2.fastq

"""

from Bio import SeqIO
import pandas as pd
import sys
import os
from collections import Counter
import gzip
from mimetypes import guess_type
from functools import partial

umi_cutoff = int(sys.argv[1])
umi_table  = sys.argv[2]
in_fastq1  = sys.argv[3]
in_fastq2  = sys.argv[4]
out_fastq1 = sys.argv[5]
out_fastq2 = sys.argv[6]

# read in reads per umi df:
reads_per_umi = pd.read_csv(umi_table, header = None)
reads_per_umi.columns = ["umi", "reads"]

# calculate the average read quality for each read (R1 + R2), save to dict avg_quality:
encoding = guess_type(in_fastq1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
avg_quality = {}
with _open(in_fastq1) as f1, _open(in_fastq2) as f2:
    for record1, record2 in zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')):
        umi = record1.name.split("_")[1]
        if reads_per_umi[reads_per_umi["umi"] == umi]["reads"].iloc[0] >= umi_cutoff:
            # print(record)
            sum_1 = sum(record1.letter_annotations["phred_quality"])
            sum_2 = sum(record2.letter_annotations["phred_quality"])
            avg_qual = (sum_1 + sum_2) / (len(record1) + len(record2))
            if not umi in avg_quality:
                avg_quality[umi] = [avg_qual, record1, record2]
            else:
                if avg_qual > avg_quality[umi][0]:
                    avg_quality[umi] = [avg_qual, record1, record2]

# save output
with open(out_fastq1, "w") as f1, open(out_fastq2, "w") as f2:
    for key in avg_quality:
        SeqIO.write(avg_quality[key][1], f1, "fastq")
        SeqIO.write(avg_quality[key][2], f2, "fastq")