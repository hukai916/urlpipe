#!/usr/bin/env python

"""
Retrieve seq given ref fasta, direction, number of nt to extract, and wether to rc.

Usage:

python get_seq.py ref.fasta start|end 20 rc (for rc or no)

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fasta = sys.argv[1]
direction = sys.argv[2]
num = int(sys.argv[3])
rc = sys.argv[4]

encoding = guess_type(fasta)[1]

_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

# output snp1, snp2, and undetermined reads
# undetermined: if both snp1 and snp2, or if neither snp1 nor snp2
with _open(fasta) as f_in:
    for record in SeqIO.parse(f_in, "fasta"):
        if direction == "start":
            seq = record.seq[:num]
        elif direction == "end":
            seq = record.seq[-num:]
        else:
            exit("Not valid direction")

        if rc == "rc":
            seq = seq.reverse_complement()
        elif rc == "no":
            seq = seq
        else:
            exit("Not valid rc option")
        print(seq)
        exit()