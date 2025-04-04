# Plasmid ID

For identifying plasmid IDs and deconvoluting PCR1 barcodes in plasmid ID sequencing protocol.

## Setup

### Install `miniforge3`

The easiest way to install the necessary dependencies for this workflow is to use [Miniforge3](https://conda-forge.org/miniforge/). Download and installation instructions can be found [here](https://github.com/conda-forge/miniforge#download).

### Create a `conda` environment and install packages

Assuming you have installed and configured `miniforge3` correctly, you can create an environment with the necessary dependencies using a single command:

```bash
conda env create -n plasmid -c conda-forge -c bioconda vsearch numpy pandas biopython
```

Or you can create the environment using the `environment.yml` file included in this repository:

```bash
conda env create -n plasmid -f environment.yml
```

## Running the workflow

### Activate the `conda` environment

```bash
conda activate plasmid
```

### Merging paired end reads with `vsearch`

The workflow assumes you have paired-end reads for each sequencing sample that you wish to deconvolute. Merging of paired-end reads is done with `vsearch` via the `merge.sh` script. Running the script without any options will produce the help menu.

```bash
./merge.sh
```

```
USAGE: ./merge.sh R1.fastq.gz R2.fastq.gz
```

You will have to run the script once for each pair of reads (i.e. once per sample) and the output will be a single `<sample_name>.merged.fastq.gz` file ready for processing with the Python script.

### 1. Identifying random barcodes in a new plasmid library and creating a barcode dictionary

To create a barcode dictionary for a new plasmid ID library, use the `plasmid_make_db.py` script. Running the script with the `-h` flag will show the help menu.

```bash
./plasmid_make_db.py -h
```

```
usage: plasmid_make_db.py [-h] [-f FW_PRIMERS] [-r RV_PRIMERS] [-5 LEFT] [-3 RIGHT] [-m MIN_COUNT] [-p MIN_PURITY]
                          seqfile

Identify random plasmid ID barcodes in sequence reads and create a dictionary.

positional arguments:
  seqfile               reads.fastq.gz

options:
  -h, --help            show this help message and exit
  -f FW_PRIMERS, --fw-primers FW_PRIMERS
                        FASTA file with forward primers (default: FW_primers.fa)
  -r RV_PRIMERS, --rv-primers RV_PRIMERS
                        FASTA file with reverse primers (default: RV_primers.fa)
  -5 LEFT, --left LEFT  5-prime (left) flanking sequnce (default: TGAACTGTACAAATGAAGGT)
  -3 RIGHT, --right RIGHT
                        3-prime (right) flanking sequence (GCTT + N12 experiment tag) (default: GCTTTGTATCTTCACC)
  -m MIN_COUNT, --min-count MIN_COUNT
                        minimum read count per well (default: 0)
  -p MIN_PURITY, --min-purity MIN_PURITY
                        minimum relative abundance for a barcode in a well (default: 0.5)
```

If successful, the script will create a table of all ASVs identified in each well (`*_asv_table.csv`) and a barcode dictionary file resulting from filtering and naming the ASVs (`*_db.csv`), as well as some text output. For example:

```bash
./plasmid_make_db.py S216.merged.fastq.gz
```

```
Processed 385483 reads from S216.merged.fastq.gz
310394 reads (80.52%) matched expected read architecture
Wrote counts for 1584 unique barcodes to S216_asv_table.csv
Wrote 67 barcodes to S216_db.csv
```

Barcodes in the fincal dictionary will be named according to the well with the highest count of that barcode.

### 2. Identifying and counting barcodes in a sample based on an existing dictionary

The `plasmid_ID.py` script takes a single argument - the reads to analyze – and produces an output CSV file with the count of each barcode identified. Running the script with the `-h` flag will show the help menu.

```bash
./plasmid_ID.py -h
```

```
usage: plasmid_ID.py [-h] [-f FW_PRIMERS] [-r RV_PRIMERS] [-5 LEFT] [-3 RIGHT] seqfile dictionary

Identify plasmid ID barcodes in sequence reads and search against barcode dictionary

positional arguments:
  seqfile               reads.fastq.gz
  dictionary            barcode_dict.csv

options:
  -h, --help            show this help message and exit
  -f FW_PRIMERS, --fw-primers FW_PRIMERS
                        FASTA file with forward primers (default: FW_primers.fa)
  -r RV_PRIMERS, --rv-primers RV_PRIMERS
                        FASTA file with reverse primers (default: RV_primers.fa)
  -5 LEFT, --left LEFT  5-prime (left) flanking sequnce (default: TGAACTGTACAAATGAAGGT)
  -3 RIGHT, --right RIGHT
                        3-prime (right) flanking sequence (GCTT + N12 experiment tag) (default: GCTTTGTATCTTCACC)
```

If successful, the script will create a results file with the count for each well of each barcode matched in the dictionary, and will produce some messages as output. For example, using the dictionary created in the previous step:

```bash
./plasmid_ID.py S216.merged.fastq.gz S216_db.csv
```

```
Processed 385483 reads from S216.merged.fastq.gz
310394 reads (80.52%) matched expected read architecture
Wrote counts for 68 matched barcodes to S216_results.csv
```

### License information

This repository contains code from [marcelamendoza/Plasmid-ID](https://github.com/marcelamendoza/Plasmid-ID) in accordance with the [MIT license](https://github.com/marcelamendoza/Plasmid-ID?tab=MIT-1-ov-file#readme), Copyright (c) 2018 marcelamendoza.

