[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/komb/README.html) [![Anaconda-Server Badge](https://anaconda.org/bioconda/komb/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/komb/badges/downloads.svg)](https://anaconda.org/bioconda/komb) [![Anaconda-Server Badge](https://anaconda.org/bioconda/komb/badges/platforms.svg)](https://anaconda.org/bioconda/komb) [![Anaconda-Server Badge](https://anaconda.org/bioconda/komb/badges/latest_release_date.svg)](https://anaconda.org/bioconda/komb)


![KOMB Logo](Images/Images_Logo.png)
## KOMB
#### Welcome to the KOMB repo! Get ready to KOMB through your (large) metagenomes to find regions of biological (functional or taxonomic) interest!
#### KOMB version: 2.0.0

## KOMB installation

Current version of KOMB has to be installed from source. It has only been tested on Linux systems. 

In order to install KOMB you will need to install several dependencies first. We recommend using conda for managing KOMB dependencies. Below is an example of installing required tools with conda.

```

```

## Quickstart example

You can test KOMB by running the following command:    
`komb -r example/2bact_42.read1.fq,example/2bact_42.read2.fq -k 51 -l 100 -t 10`    

## KOMB usage
Installation will create  an executable `komb` inside it. Once the binary is obtained you can add it to path for ease of use. Users can run KOMB as follows using various command line options.

```
USAGE: 

   komb  [-o <string>] [-t <int>] [-k <int>] [-n <int>] [-l <int>] -r
               <string> [-a] [-f] [-b] [-s] [--] [--version] [-h]


Where: 

   -o <string>,  --output <string>
     Output directory [Default: output_yyyymmdd_hhmmss]

   -t <int>,  --threads <int>
     Number of Threads [Default: Max]

   -k <int>,  --kmer <int>
     Kmer size for Abyss, Bifrost uses 31 [Default: 31]

   -n <int>,  --numhits <int>
     Bowtie2 maximum number of hits [Default: 1000]

   -l <int>,  --readlen <int>
     Read Length (can be average) [Default: 100]

   -r <string>,  --reads <string>
     (required)  Paired-read file separated by ',' [Default: read1.fq
     ,read2.fq]

   -a,  --alignment
     Keep alignment files [Default: delete alignment]

   -f,  --fasta
     Reads provided are fasta files [Default: fastq]

   -b,  --bifrost
     Run bifrost instead of abyss [Default: run abyss]

   -s,  --spades
     Runs spades and uses GFA graph instead of bifrost + bowtie2 [Default:
     run abyss]

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   KOMB: Taxonomy-oblivious characterization of metagenome dynamics
   ```
