# Plasmid ID

For identifying plasmid IDs and deconvoluting PCR1 barcodes in plasmid ID sequencing protocol.

## Setup

### Install `vsearch`

First install `vsearch` according to the [install documentation](https://github.com/torognes/vsearch#download-and-install) for your particular operating system.

### Install Python dependencies

To manage Python dependencies, it is recommended to use Python virtual environments, which is the method used below. `conda` or other Python package managers may also be used if preferred.

The following commands create a virtual environment, activate the environment, and install the necessary Python packages.

```bash
python -m venv plasmid_env
source plasmid_env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

You are now ready to run the workflow.

## Running the workflow

### Merging paired end reads with `vsearch`

The workflow assumes you have paired-end reads for each sequencing sample that you wish to deconvolute. This step is accomplished using the `merge.sh` script. Running the script without any options will produce the help menu.

```bash
./merge.sh
```

```
USAGE: ./merge.sh R1.fastq.gz R2.fastq.gz
```

You will have to run the script once for each pair of reads (i.e. once per sample) and the output will be a single `<sample_name>.merged.fastq.gz` file ready for processing with the Python script.

### Identifying and counting barcodes








