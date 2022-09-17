"""
To annotate reads into no_indel, indel_5p, indel_3p, and indel_5p_3p (flanking)
    # align R1 to 20bp before CAG repeat, allow 1 mismatch: TCGAGTCCCTCAAGTCCTTC
    # align R2 to 20bp after CAG repeat, allow 1 mismatch: CCGCCACCGCCGCCGCCGCC (use rev comp)
    # if both match: "no indel"; if R1 doesn't match: "indel, 5p"; if R2 doesn't match: "indel, 3p"
    # implementation: instead of perform mapping, simple search sub_string with python regex:
    regex.search("(xxgyy){e<=1}", "xxggyy") # https://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string, note, regex.search("xxgyy{e<=1}", "xxggyy") does not work.
Usage:
    python annotate_indel.py r1.fastq r2.fastq dir_to_store_no_indel dir_to_store_indel_5p dir_to_store_indel_3p dir_to_store_indel_5p_3p
"""

from Bio import SeqIO
import sys
import os
from collections import Counter
import regex

r1 = sys.argv[1]
r2 = sys.argv[2]
no_indel_dir = sys.argv[3]
indel_5p_dir = sys.argv[4]
indel_3p_dir = sys.argv[5]
indel_5p_3p_dir = sys.argv[6]

r1_20nt = "TCGAGTCCCTCAAGTCCTTC"
r2_20nt = "GGCGGCGGCGGCGGTGGCGG" # rev complemented of CCGCCACCGCCGCCGCCGCC
mismatch = 2 # including insertion and deletion

r1_match = {}
r2_match = {}
r_match  = {}

with open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r1_20nt + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)):
            r1_match[record.name] = 1
        else:
            r1_match[record.name] = 0
with open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r2_20nt + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)):
            r2_match[record.name] = 1
        else:
            r2_match[record.name] = 0

for k in r1_match:
    assert k in r2_match, "not equal records!"
    r_match[k] = r1_match[k] + r2_match[k]

# output:
out_no_indel_r1 = os.path.join(no_indel_dir, os.path.basename(r1))
out_no_indel_r2 = os.path.join(no_indel_dir, os.path.basename(r2))
out_indel_5p_r1 = os.path.join(indel_5p_dir, os.path.basename(r1))
out_indel_5p_r2 = os.path.join(indel_5p_dir, os.path.basename(r2))
out_indel_3p_r1 = os.path.join(indel_3p_dir, os.path.basename(r1))
out_indel_3p_r2 = os.path.join(indel_3p_dir, os.path.basename(r2))
out_indel_5p_3p_r1 = os.path.join(indel_5p_3p_dir, os.path.basename(r1))
out_indel_5p_3p_r2 = os.path.join(indel_5p_3p_dir, os.path.basename(r2))

os.makedirs(os.path.dirname(out_no_indel_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_no_indel_r2), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_5p_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_5p_r2), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_3p_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_3p_r2), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_5p_3p_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_indel_5p_3p_r2), exist_ok=True)

count_5p = 0
count_3p = 0

r1_no_indel = []
r1_indel_5p = []
r1_indel_3p = []
r1_indel_5p_3p = []
with open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r1_no_indel.append(record)
        elif r_match[record.name] == 1:
            if r1_match[record.name]:
                r1_indel_3p.append(record)
                count_3p += 1
            else:
                r1_indel_5p.append(record)
                count_5p += 1
        elif r_match[record.name] == 0:
            r1_indel_5p_3p.append(record)
SeqIO.write(r1_no_indel, out_no_indel_r1, "fastq")
SeqIO.write(r1_indel_5p, out_indel_5p_r1, "fastq")
SeqIO.write(r1_indel_3p, out_indel_3p_r1, "fastq")
SeqIO.write(r1_indel_5p_3p, out_indel_5p_3p_r1, "fastq")

r2_no_indel = []
r2_indel_5p = []
r2_indel_3p = []
r2_indel_5p_3p = []

with open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r2_no_indel.append(record)
        elif r_match[record.name] == 1:
            if r2_match[record.name]:
                r2_indel_5p.append(record)
            else:
                r2_indel_3p.append(record)
        elif r_match[record.name] == 0:
            r2_indel_5p_3p.append(record)
SeqIO.write(r2_no_indel, out_no_indel_r2, "fastq")
SeqIO.write(r2_indel_5p, out_indel_5p_r2, "fastq")
SeqIO.write(r2_indel_3p, out_indel_3p_r2, "fastq")
SeqIO.write(r2_indel_5p_3p, out_indel_5p_3p_r2, "fastq")


# print some stats:
res = Counter(r_match.values())
with open("annotation_indel_stats_percentage.tsv", "a") as f:
    name = os.path.basename(r1).split("_R1_")[0]
    p2 = str(res[2]/(sum([res[2], res[1], res[0]])))
    p_5p = str(count_5p/(sum([res[2], res[1], res[0]])))
    p_3p = str(count_3p/(sum([res[2], res[1], res[0]])))
    p0 = str(res[0]/(sum([res[2], res[1], res[0]])))

    f.write(name + "\t" + str(res[2]) + "\t" + p2 + "\t" + str(count_5p) + "\t" + p_5p + "\t" + str(count_3p) + "\t" + p_3p + "\t" + str(res[0]) + "\t" + p0 + "\n")

# f_out_htt_r1.close()
# f_out_htt_r2.close()
# f_out_mis_r1.close()
# f_out_mis_r2.close()
# f_out_problem_r1.close()
# f_out_problem_r2.close()
