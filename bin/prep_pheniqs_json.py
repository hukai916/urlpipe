#!/usr/bin/env python

"""
Prepare json for pheniqs given barcodes info.

Usage:
python prep_pheniqs_json.py sequence.fastq.gz outfile_prefix_name_prefix rc(boolean) barcode_file

"""

import json
import sys
from Bio.Seq import Seq 

seq = sys.argv[1]
outfile_prefix = sys.argv[2]
rc = sys.argv[3]
barcode_pos = sys.argv[4]
barcode_file = sys.argv[5]

if rc == 1:
    outfile_prefix += "_rc"

bc_condition = {}
with open(barcode_file) as f:
    for line in f:
        if ":" in line:
            sample_name, sample_barcode = line.strip().split(":")
            if not sample_name in bc_condition:
                sample_barcode = sample_barcode.strip()
                if rc:
                   sample_barcode = str(Seq(sample_barcode).reverse_complement())
                bc_condition[sample_name] = sample_barcode.strip()
                    
config = {}
config['report url'] = outfile_prefix + "_pheniqs_report.json"
config['input'] = [seq]
config['template'] = {}
config['template']['transform'] = {}
config['template']['transform']['token'] = ["0::"]
config['sample'] = {}
config['sample']['algorithm'] = 'pamld'
# config['sample']['confidence threshold'] = 0.99
# config['sample']['noise'] = 0.01
config['sample']['transform'] = {}
# config['sample']['transform']['token'] = ["0::24"]
config['sample']['transform']['token'] = [barcode_pos]
config['sample']['codec'] = {}

for bc_cond in bc_condition:
    name = bc_cond
    config['sample']['codec']['@' + name] = {}
    config['sample']['codec']['@' + name]['LB'] = name
    config['sample']['codec']['@' + name]['barcode'] = [bc_condition[bc_cond]]
    config['sample']['codec']['@' + name]['output']  = [outfile_prefix + "_" + name + ".fastq.gz"]
    config['undetermined'] = {}
    config['undetermined']["output"] = [outfile_prefix + "_undetermined.fastq.gz"]

with open(outfile_prefix + ".json", 'w') as outfile_prefix:
    json.dump(config, outfile_prefix, indent = 4)
