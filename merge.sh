#!/bin/bash

if [ $# -ne 2 ]; then
	echo "USAGE: $0 R1.fastq.gz R2.fastq.gz"
	exit 1
fi

R1=$1
R2=$2

printf "\nMerging $R1 and $R2\n"

# Extract sample names and verify that they match
S1=$(echo $R1 | sed -e 's/_R1_.*fastq.gz//')
S2=$(echo $R2 | sed -e 's/_R2_.*fastq.gz//')

if [ ! $S1 == $S2 ]; then
	echo "Sample names do not match!"
	exit 1
fi

# Merge reads
vsearch --fastq_mergepairs $R1 \
		--reverse $R2 \
		--fastqout ${S1}.merged.fastq \
		--fastq_allowmergestagger
		#--threads $NCPUS

gzip ${S1}.merged.fastq

