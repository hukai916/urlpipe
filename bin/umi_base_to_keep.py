#!/usr/bin/env python

"""
To trim UMI that is appended to the read ID.

Usage:
python umi_base_to_keep.py 12 output.fastq.gz

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]
umi_base_to_keep = sys.argv[2]
fastq_out = sys.argv[3]

encoding = guess_type(fastq)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
_open_out = partial(gzip.open, mode='wt') if encoding == 'gzip' else open


with _open(fastq) as f_in, _open_out(fastq_out) as f_out:
    for record in SeqIO.parse(f_in, "fastq"):
        print(dir(record))
        print(record.id)
        print(record.name)
        exit()
