#!/usr/bin/env python

import sys 
from Bio import SeqIO

ref_fa = sys.argv[1]
start, end = int(sys.argv[2]), int(sys.argv[3])
repeat_unit = sys.argv[4]
repeat_range = sys.argv[5]
out_ref = sys.argv[6]

fasta_sequences = SeqIO.parse(ref_fa, "fasta")

for fasta in fasta_sequences:
    id = fasta.id
    seq = fasta.seq
    before_repeat = seq[:start - 1]
    after_repeat = seq[end:]

repeat_range_s, repeat_range_e = repeat_range.split(":")[:2]

fa = []
for i in range(int(repeat_range_s), int(repeat_range_e) + 1):
    id = "ref_repeat_" + str(i * 3)
    seq = before_repeat + repeat_unit.upper() * i + after_repeat
    fa.append(SeqIO.SeqRecord(seq = seq, id = id, description = "ref_with_" + str(i) + "_repeat_units"))

with open(out_ref, "w") as output_handle:
    SeqIO.write(fa, output_handle, "fasta")