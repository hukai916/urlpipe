#!/usr/bin/env python

"""
Plot UMI distribution.
"""

import sys
import os
from collections import Counter
from statistics import mean
from utils import plot_repeat_dist

tsv         = sys.argv[1]
sample_name = sys.argv[2]
outdir      = sys.argv[3]
cutoff      = sys.argv[4]

with open(outdir + "/test.txt", "w") as f:
    f.write("Simply a test!")
