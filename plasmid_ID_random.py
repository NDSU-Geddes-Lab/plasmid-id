#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import csv
import gzip
import itertools
import pandas as pd
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser(
                    prog='plasmid_ID_random.py',
                    description='Identify random plasmid ID barcodes in sequence reads')

parser.add_argument('seqfile', help='reads.fastq.gz')
parser.add_argument('-f', '--fw-primers',
                    help='FASTA file with forward primers',
                    default='FW_primers.fa')
parser.add_argument('-r', '--rv-primers',
                    help='FASTA file with reverse primers',
                    default='RV_primers.fa')
parser.add_argument('-5', '--left',
                    help='5-prime (left) flanking sequnce',
                    default='TGAACTGTACAAATGAAGGT')
parser.add_argument('-3', '--right',
                    help='3-prime (right) flanking sequence (GCTT + N12 experiment tag)',
                    default='GCTTTGTATCTTCACC')

args = parser.parse_args()

# Reead in forward and reverse primer sequences
forward_dict = SeqIO.to_dict(SeqIO.parse(args.fw_primers, "fasta"))
forward_dict = {k:str(v.seq) for k, v in forward_dict.items()}
reverse_dict = SeqIO.to_dict(SeqIO.parse(args.rv_primers, "fasta"))
reverse_dict = {k:str(v.seq.reverse_complement()) for k, v in reverse_dict.items()}

# Create a dictionary to store each identified barcode
barcodes = {}
for col in reverse_dict.keys():
    for row in forward_dict.keys():
        well = col+row
        barcodes[well] = {}

# Identify well and barcode all at once
def find_barcode_and_well(seq, fwd_primers, rev_primers):
    """Iterate through forward and reverse primer pairs
    and identify well and barcode.
    """
    for row, fwd_seq in fwd_primers.items():
        for col, rev_seq in rev_primers.items():
            # Read architecture: {fwd}{right}{N33}{left}{rev}
            pattern = fwd_seq + args.left + "(.{33})" + args.right + rev_seq
            hit = re.search(pattern, seq)
            if hit is not None:
                return(col+row, hit.group(1))
    return(None, None)

# Main loop to iterate through each fastq sequence
with gzip.open(args.seqfile, "rt") as r1:
    for fw in SeqIO.parse(r1, "fastq") :
        str_seq = str(fw.seq)
        well, barcode = find_barcode_and_well(str_seq, forward_dict, reverse_dict)
        # If well is None, that means no barcode was found
        if well is None:
            continue
        # If we haven't seen this barcode in this well, then we set count to 1
        elif barcode not in barcodes[well]:
            barcodes[well][barcode] = 1
        # If we've seen it already in that well, then increment the count
        else:
            barcodes[well][barcode] += 1

# Convert to dataframe and output as CSV
sample_name = args.seqfile.split('.')[0]
outfile = sample_name + "_barcodes.csv"
results = pd.DataFrame.from_dict(barcodes).fillna(0).astype('int')
n_barcodes = len(results)
#print(results)
results.to_csv(outfile)
print(f"Wrote counts for {n_barcodes} unique barcodes to {outfile}")

