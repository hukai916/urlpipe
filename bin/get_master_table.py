#!/usr/bin/env python

import sys
import pandas as pd

frac_csv = sys.argv[1]
indel_csv = sys.argv[2]
allele_number = sys.argv[3]
outfile = sys.artv[4]

df_frac = pd.read_csv(frac_csv)
df_indel = pd.read_csv(indel_csv)

with open(outfile, "w") as f:
    f.write("test")
    # df_res.to_csv("name.csv", index = False)