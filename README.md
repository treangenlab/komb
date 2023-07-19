![KOMB Logo](Images/Images_Logo.png)
## KOMB
#### Welcome to the KOMB repo! Get ready to KOMB through your (large) metagenomes to find regions of biological (functional or taxonomic) interest!
#### KOMB version: 2.0.0

## KOMB installation

Current version of KOMB has to be installed from source. It has only been tested on Linux systems. 

In order to install KOMB you will need to install several dependencies first. We recommend using `conda` ([Miniconda download](https://docs.conda.io/en/latest/miniconda.html)) for managing KOMB dependencies. Below is an example of installing required tools with `conda`.

1. Create a new `conda` environment and make sure that `conda-forge` and `bioconda` channels are enabled.
```
conda create --name komb-env python=3.9
conda config --add channels conda-forge
conda config --add channels bioconda
```
2. Install dependencies available through `conda`
```
conda install bwa-mem2
conda install seqkit
conda install igraph>=0.10.0
```
**Note:** we have upgraded KOMB to be compatible with newer version of the [igraph library](https://igraph.org/c/) which means that the versions 0.8.3 and older no longer will work. We hope that igraph API will stay consitent going forward, but we have no way to ensure that.
3. Install `rustup` and `nightly` toolchain and build `ggcat` from source
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup toolchain install nightly
git clone https://github.com/algbio/ggcat --recursive
cd ggcat/
cargo install --path crates/cmdline/ --locked
```
4. Clone the repo and build the source code
```
git clone https://gitlab.com/treangenlab/komb.git.
cd komb 
./autogen.sh
./configure
make; make install
```
5. Now you will have a `komb2` executable in the `komb/bin` directory and you should be able to run `KOMB.py` in order to run the complete KOMB pipeline.

## Quickstart example

You can test KOMB by running the following command:    
`KOMB.py -i example_data/reads1.fastq -j example_data/reads2.fastq -k 51`

Example data is hosted on OSF at the following URL: [ADD URL] and on Zenodo under DOI: [ADD DOI].

## KOMB usage
Full set of parameters available in the KOMB pipeline is shown below and can be accessed by running `python KOMB.py --help`.

```
usage: KOMB.py [-h] -i INPUT_READS1 -j INPUT_READS2 [-o OUTPUT_DIR] [--keep-alignments] [-e LOG_FILE] [--overwrite] -k KMER_SIZE [-t NUM_THREADS] [-l MIN_UNITIG_LENGTH] [-v VERBOSITY] [-c MIN_COUNT]
               [-m GGCAT_MEMORY] [--eulertigs | --greedy-matchtigs | --pathtigs] [--min-seed-length MIN_SEED_LENGTH]

KOMB Analysis Pipeline example: python KOMB.py -i <read1.fq> -j <read2.fq> -k <k-mer size> -t <threads>

optional arguments:
  -h, --help            show this help message and exit
  --eulertigs           Generate Eulertigs instead of unitigs
  --greedy-matchtigs    Generate greedy matchtigs instead of unitigs
  --pathtigs            Generate pathtigs instead of unitigs

Input/Output:
  -i INPUT_READS1, --input-reads1 INPUT_READS1
                        Path to the first sequencing reads file for paired-end data in FASTQ format
  -j INPUT_READS2, --input-reads2 INPUT_READS2
                        Path to the second sequencing reads file for paired-end data in FASTQ format
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Path to the first sequencing reads file for paired-end data in FASTQ format
  --keep-alignments     Keep SAM files after the graph has been constructed (might require a lot of disk space)
  -e LOG_FILE, --log-file LOG_FILE
                        File for logging [default: stdout]
  --overwrite           Delete <output-dir> and create a new one in case it exists

Common arguments:
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        k-mer size used for the *tig construction and subsequent analyses
                        use -1 to let KOMB automatically pick a value [default: -1]
  -t NUM_THREADS, --num-threads NUM_THREADS
                        Maximum number of threads you want programs to use, note that some might use less than the amount specified
  -l MIN_UNITIG_LENGTH, --min-unitig-length MIN_UNITIG_LENGTH
                        Minimum length of a unitig to be kept for the analysis.
                        Value -1 indicates setting this to the read length [default]
                        Value 0 would result in keeping all unitigs, and values > 0 will apply the filter
  -v VERBOSITY, --verbosity VERBOSITY
                        Logging level: 0 (DEBUG), 1 (INFO), 2 (ERROR)

GGCAT *tig construction:
  -c MIN_COUNT, --min-count MIN_COUNT
                        Minimum count required to keep a kmer [default: 2]
  -m GGCAT_MEMORY, --ggcat-memory GGCAT_MEMORY
                        Maximum memory usage for GGCAT (GB) [default: 8]

BWA MEM parameters:
  --min-seed-length MIN_SEED_LENGTH
                        Minimum seed length. Matches shorter than the value will be missed.
```
