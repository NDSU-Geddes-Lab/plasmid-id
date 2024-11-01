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

### Identifying and counting barcodes with `plasmid_ID.py`

The `plasmid_ID.py` script takes a single argument - the reads to analyze – and produces an output CSV file with the count of each barcode identified. Running the script with the `-h` flag will show the help menu.

```bash
./plasmid_ID.py -h
```

```
usage: plasmid_ID.py [-h] seqfile

Identifies plasmid ID barcodes in sequence reads

positional arguments:
  seqfile     reads.fastq.gz

options:
  -h, --help  show this help message and exit
```

For example, to process a merged sequence file from the previous step:

```bash
./plasmid_ID.py sample1.merged.fastq.gz
```

If successful, the script will not produce any text output, but you should see two new files created:

```
sample1_distances_merge.txt
sample1_results.csv
```

