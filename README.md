![alt text](Images/Logo.png)
# KOMB
KOMB V1.0

### Dependencies
* Bifrost: `conda install -c bioconda bifrost`
* Abyss: `conda install -c bioconda abyss`
* igraph: `conda install -c conda-forge igraph`
* Bowtie2: `conda install -c bioconda bowtie2`


### Installation
KOMB is installed through autotools. We've tested KOMB on both Linux and OSX.

##### On OSX:

`./autogen.sh`   
`./configure CXX="/usr/local/Cellar/llvm/10.0.1_1/bin/clang++"`   
`make`   
`make install`

Please ensure your CXX supports openmp (-fopenmp). You can install the latest llvm through brew: `brew install llvm` that should be compatible.

**Important:** Users must also check if the path to the miniconda installation in `src/Makefile.am` is accurate and make any changes as needed.

##### On Linux:

`./autogen.sh`   
`./configure`   
`make`   
`make install`

### KOMB usage
Installation will create a `bin` folder and an executable `komb` inside it. Users can run KOMB as follows using various command line options.

```
USAGE: 

   ./bin/komb  [-o <string>] [-t <int>] [-k <int>] [-n <int>] [-l <int>] -r
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
