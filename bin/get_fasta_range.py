#!/usr/bin/env python

"""
Return base in given range.

Usage:
python get_fasta_range.py sample.fasta 10:40

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fasta = sys.argv[1]
range_to_return = sys.argv[2]

encoding = guess_type(fasta)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

tem_range = [int(x) if x != '' else x for x in range_to_return.split(":")]

with _open(fasta) as f:
    for record in SeqIO.parse(f, "fasta"):
        if tem_range[0] == '' and tem_range[1] == '':
            print(record.seq.upper())
        elif tem_range[0] == '':
            print(record.seq.upper()[:tem_range[1]])
        elif tem_range[1] == '':
            print(record.seq.upper()[tem_range[0]:])
        else:
            print(record.seq.upper()[tem_range[0]:tem_range[1]])
