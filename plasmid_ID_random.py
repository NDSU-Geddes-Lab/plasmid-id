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
import sys

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
parser.add_argument('-m', '--min-count',
                    help='minimum read count per well',
                    type=int, default=0)
parser.add_argument('-p', '--min-purity',
                    help='minimum relative abundance for a barcode in a well',
                    type=float, default=0.5) # Ensures at most 1 barcode will be kept per well

args = parser.parse_args()
min_count = args.min_count
min_purity = args.min_purity
if min_purity < 0.5:
    print("WARNING: setting --min-purity < 0.5 may result in multiple barcodes per well in your dictionary!")

# Read in forward and reverse primer sequences
forward_dict = SeqIO.to_dict(SeqIO.parse(args.fw_primers, "fasta"))
forward_dict = {k:str(v.seq) for k, v in forward_dict.items()}
reverse_dict = SeqIO.to_dict(SeqIO.parse(args.rv_primers, "fasta"))
reverse_dict = {k:str(v.seq.reverse_complement()) for k, v in reverse_dict.items()}

# Create a dictionary to store each identified barcode
plate = {}
for col in reverse_dict.keys():
    for row in forward_dict.keys():
        well = col+row
        plate[well] = {}

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
n_reads = 0
n_matched = 0
with gzip.open(args.seqfile, "rt") as r1:
    for fw in SeqIO.parse(r1, "fastq") :
        n_reads += 1
        str_seq = str(fw.seq)
        well, barcode = find_barcode_and_well(str_seq, forward_dict, reverse_dict)
        # If well is None, that means no barcode was found
        if well is None:
            continue
        # If we haven't seen this barcode in this well, then we set count to 1
        elif barcode not in plate[well]:
            plate[well][barcode] = 1
            n_matched += 1
        # If we've seen it already in that well, then increment the count
        else:
            plate[well][barcode] += 1
            n_matched += 1

# Print some basic diagnostics
pct_matched = round((n_matched/n_reads)*100, 2)
print(f"Processed {n_reads} reads from {args.seqfile}")
print(f"{n_matched} reads ({pct_matched}%) matched expected read architecture")

# Output entire ASV table for reference, before we start filtering anything out
# Convert to dataframe and output as CSV
sample_name = args.seqfile.split('.')[0]
outfile = sample_name + ".asv_table.csv"
results = pd.DataFrame.from_dict(plate).fillna(0).astype('int')

# Create output and report results
n_barcodes = len(results)
results.to_csv(outfile)
print(f"Wrote counts for {n_barcodes} unique barcodes to {outfile}")

# Now we can proceed with creating the database
# Remove wells with less than --min-count reads
low_count = []
if min_count > 0:
    for well, bc_dict in plate.items():
        well_total = sum(bc_dict.values())
        if well_total < min_count:
            low_count.append(well)
    for well in low_count:
        plate.pop(well)
    print(f"Removed {len(low_count)} wells with less than {min_count} reads: {low_count}")
    print(f"{len(plate.keys())} wells remaining")

# Remove low purity barcodes before final deduplication
if min_purity > 0:
    for well, bc_dict in plate.items():
        low_purity = []
        well_total = sum(bc_dict.values())
        for bc, count in bc_dict.items():
            if (count/well_total) < min_purity:
                low_purity.append(bc)
        for bc in low_purity:
            bc_dict.pop(bc)

# Now we're done with per-well analysis, so we need to rekey on barcode sequence
barcodes = {}
for well, bc_dict in plate.items():
    for bc, count in bc_dict.items():
        if bc not in barcodes:
            barcodes[bc] = {}
        barcodes[bc][well] = count

# Now iterate through barcodes and write the top count well for each to the db.csv
db_file = sample_name + "_db.csv"
db = open(db_file, "w")
db.write("barcode,sequence\n")
n_barcodes_final = 0
for bc_seq, well_dict in barcodes.items():
    top_well = max(well_dict, key=well_dict.get)
    bc_name = sample_name + "_" + top_well.lower()
    db.write(f"{bc_name},{bc_seq}\n")
    n_barcodes_final += 1

db.close()
print(f"Wrote {n_barcodes_final} barcodes to {db_file}")

