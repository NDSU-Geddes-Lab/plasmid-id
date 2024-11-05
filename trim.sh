#!/bin/bash

if [ $# -eq 1 ]; then
	FASTQ=$1
	LENGTH=150
elif [ $# -eq 2 ]; then
	FASTQ=$1
	LENGTH=$2
else
	echo "USAGE: $0 reads.fastq.gz [length]"
	exit 1
fi

# Extract sample names
SAMPLE=$(echo $FASTQ | cut -d. -f1)

# Truncate reads to 150 bp
#R1
vsearch --fastx_filter $FASTQ \
		--fastqout ${SAMPLE}.trimmed.fastq \
		--fastq_trunclen $LENGTH

