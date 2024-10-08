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


# Here we read all the sequence indeces for the forward and reverse reads

# In[2]:


forward_dict = SeqIO.to_dict(SeqIO.parse("FW_primers.fa", "fasta"))
forward_dict = {k:str(v.seq) for k, v in forward_dict.items()}
reverse_dict = SeqIO.to_dict(SeqIO.parse("RV_primers.fa", "fasta"))
reverse_dict = {k:str(v.seq.reverse_complement()) for k, v in reverse_dict.items()}


# Now, we read all the plasmid sequences

# In[3]:


plasmids = {}
with open('Plasmid_ID_correct_position_jake_rcbc100.csv', 'r') as csvfile:
    table_reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in table_reader:
        plasmids[row["Barcode name"]]= row["ID plasmid sequence"]


# In[4]:


table_reader


# Let's iterate the two sequence files and produce the table with the distances

# In[6]:


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

fw_f = "../Jake13.merged.fastq.gz"
distances_f = "distances_merge.txt"
with gzip.open(fw_f, "rt") as r1,open(distances_f,"w") as f3:
    f3.write("fw_name,fw_pos,rev_name,rev_pos,pl_name,pl_pos\n")
    for fw in SeqIO.parse(r1, "fastq") :
        #print(fw)
        str_seq = str(fw.seq)
        fw_name, fw_pos = find_hash_position(str_seq, forward_dict)
        rv_name, rv_pos = find_hash_position(str_seq, reverse_dict)
        pl_name, pl_pos = find_hash_position(str_seq, plasmids)
        f3.write(",".join([fw_name, str(fw_pos), rv_name, str(rv_pos), pl_name,str(pl_pos)])) 
        f3.write("\n")
    
    


# In[7]:


df = pd.read_csv(distances_f)
df


# In[8]:


summ=df[['fw_name', 'rev_name', 'pl_name', 'pl_pos']].groupby(['fw_name', 'rev_name', 'pl_name'], as_index=False).agg( ['count','mean'])

summ.columns = ['fw_name','rev_name','pl_name','count','mean']

# In[9]:


summ.to_csv("results.csv", sep=',')


# In[10]:


summ


# In[ ]:




