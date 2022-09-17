"""
To create an artifical reference genome: CAG length ranging from 0 - 300nt.

Usage:

python create_reference.py outfile
"""

import sys
# output_dir = sys.argv[1]

# output_file = os.path.join(output_dir, "ref_cag_500.fasta")
# os.makedirs(os.path.dirname(output_file), exist_ok=True)

part1  = "CCCATCGGGCAGGAAGCCGTCATGGCAACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC"
part2 = "CCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGC"
cag_450 = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
# original is 150, the second last is caa, extend 2 times while fixing caa to cag for the extended nt.


for i in range(1, 451):
    print(">CAG_" + str(i))
    print(part1 + cag_450[0: i] + part2)
