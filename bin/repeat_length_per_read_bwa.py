#!/usr/bin/env python

"""
Count the repeat length by parsing the BAM file. The reference id contains the length of repeat length.

"""

from Bio import SeqIO
import sys
import pysam
from mimetypes import guess_type
from functools import partial
import csv
import gzip
import numpy as np

bam = sys.argv[1]
reads = sys.argv[2]
outfile = sys.argv[3]

encoding = guess_type(reads)[1]  # uses file extension
_open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open

reads_dict = {}

with _open(reads) as f:
    for record in SeqIO.parse(f, 'fastq'):
       reads_dict[record.id] = np.nan

with pysam.AlignmentFile(bam, "rb") as bam:
    for read in bam:
        if not read.reference_name == None:
            reads_dict[read.query_name] = read.reference_name.split("_")[-1]

with open(outfile, "w", newline = "") as f:
    writer = csv.writer(f)
    for key, value in reads_dict.items():
        writer.writerow([key, value, np.nan])