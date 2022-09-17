"""
To count the CAG repeat length.

Usage:
    python count_cag.py r1.fastq dir_output
"""

from Bio import SeqIO
import sys
import os
from collections import Counter
import regex

r1 = sys.argv[1]
output_dir = sys.argv[2]

r1_20nt = "TCGAGTCCCTCAAGTCCTTC"
r2_20nt = "CCGCCACCGCCGCCGCCGCC"
mismatch = 2
dict_cag = {"problem": 0, "plus": 0}

with open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        r1_search = regex.search("(" + r1_20nt + ")" + "{e<=" + str(mismatch) + "}", str(record.seq))
        r2_search = regex.search("(" + r2_20nt + ")" + "{e<=" + str(mismatch) + "}", str(record.seq))
        if not r1_search:
            dict_cag["problem"] += 1 # should not happen
        elif not r2_search:
            dict_cag["plus"] += 1
            # print(record.seq)
        else:
            r1_match_start = r1_search.start()
            r1_match_length = len(r1_search.captures()[0])
            r2_match_end = r2_search.end()
            r2_match_length = len(r2_search.captures()[0])
            cag_length = r2_match_end - r1_match_start - r1_match_length - r2_match_length
            # compensate for the insersion match of r2:
            if r2_search.captures()[0][:2] == "AG":
                cag_length += 2
            else if r2_search.captures()[0][:1] == "A":
                cag_length += 1

            if not cag_length in dict_cag:
                dict_cag[cag_length] = 1
            else:
                dict_cag[cag_length] += 1
keys = [x for x in dict_cag.keys() if not x == "problem" and not x == "plus"]

# for i in sorted(keys):
#     print(i, dict_cag[i])
# print("plus", "\t", dict_cag["plus"])
# print("problem", "\t", dict_cag["problem"])

# output:
output_file = os.path.join(output_dir, os.path.basename(r1) + ".tsv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)
print(output_file)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "\t" + str(dict_cag[i]) + "\n")
    f.write("plus" + "\t" + str(dict_cag["plus"]) + "\n")
    f.write("problem" + "\t" + str(dict_cag["problem"]) + "\n")
