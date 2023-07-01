#!/usr/bin/env python

"""
Return fastq read count.

Usage:
python get_fastq_count.py sample.fastq.gz

"""

from Bio import SeqIO
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]

encoding = guess_type(fastq)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

count = 0
with _open(fastq) as f:
    for record in SeqIO.parse(f, "fastq"):
        count += 1

print(count)