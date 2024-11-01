#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.Seq import Seq
import csv
import gzip
import itertools
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(
                    prog='plasmid_ID.py',
                    description='Identifies plasmid ID barcodes in sequence reads')

parser.add_argument('seqfile', help='reads.fastq.gz')
args = parser.parse_args()

# Here we read all the sequence indeces for the forward and reverse reads
forward_dict = SeqIO.to_dict(SeqIO.parse("FW_primers.fa", "fasta"))
forward_dict = {k:str(v.seq) for k, v in forward_dict.items()}
reverse_dict = SeqIO.to_dict(SeqIO.parse("RV_primers.fa", "fasta"))
reverse_dict = {k:str(v.seq.reverse_complement()) for k, v in reverse_dict.items()}


# Now, we read all the plasmid sequences
plasmids = {}
with open('rcbc100.csv', 'r') as csvfile:
    table_reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in table_reader:
        plasmids[row["Barcode name"]]= row["ID plasmid sequence"]

table_reader

# Iterate the two sequence files and produce the table with the distances
def find_hash_position(sequence, tags):
    #print("________")
    #print(sequence)
    #print("________")
    for key, value in tags.items():
    #    print(value)
        position = sequence.find(value)
        if position != -1:
            return (key, position)
    return ("none", -1)

sample_name = args.seqfile.split('.')[0]
distances_f = sample_name + "_distances_merge.txt"
with gzip.open(args.seqfile, "rt") as r1,open(distances_f,"w") as f3:
    f3.write("fw_name,fw_pos,rev_name,rev_pos,pl_name,pl_pos\n")
    for fw in SeqIO.parse(r1, "fastq") :
        #print(fw)
        str_seq = str(fw.seq)
        fw_name, fw_pos = find_hash_position(str_seq, forward_dict)
        rv_name, rv_pos = find_hash_position(str_seq, reverse_dict)
        pl_name, pl_pos = find_hash_position(str_seq, plasmids)
        f3.write(",".join([fw_name, str(fw_pos), rv_name, str(rv_pos), pl_name,str(pl_pos)])) 
        f3.write("\n")
    
df = pd.read_csv(distances_f)
df

summ=df[['fw_name', 'rev_name', 'pl_name', 'pl_pos']].groupby(['fw_name', 'rev_name', 'pl_name'], as_index=False).agg( ['count','mean'])

summ.columns = ['fw_name','rev_name','pl_name','count','mean']

outfile = sample_name + "_results.csv"
summ.to_csv(outfile, sep=',')

