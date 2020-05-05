![alt text](Images/Logo.png)
# KOMB - Analyzing K(c)ore repeats in (Meta)genomes
KOMB uses K-Core decomposition of unitig graphs  to analyze repeats in single and metagenomes.    

Please do go through the [Project WIKI](https://gitlab.com/treangenlab/komb/-/wikis/Getting-Started) before you use KOMB to understand the intricacies of the various parameters and modes involved.

 

## Introduction
<p style="text-align: center;">

#### Motivation
KOMB is a scalable and memory efficient tool that identify repetitive DNA in metagenomes. It is grounded on an efficient parallel unitig graph construction and uses K-Core decomposition
,a popular hierarchial graph decomposition tool to identify repetitive unitigs. Efficient identification of homoilogous regions in metagenomes has been a major challenge for the community over the past few years. Recent advances in high-throughput sequencing has cut down the time and cost of sequencing large samples resulting in the availability of a large number of metagenomic data. This is a difficult problem as metagenomic data often containes both inter- and intra-genomic sequence similarities that cause irregularities and induce errors in metagenomic sequence assembly. 

#### Previous approaches
Most current assemblers use a contig graph that collapse repetitive contig information and then use heuristics to identify highly 'tangled' nodes as repetitive. Previous approaches have included using Betweenness centrality based approaches which is computationally intensive and not scalable. Another recent approach has been using Approximate betweenness centrality which improves runtime but as an approximate method is combined with other 
features like contig length and coverage as inputs to Random Forest Model for prediction. Though this method is specific it results in sub-optimal sensitivity.

#### KOMB
In contrast, we present KOMB which uses a unitig graph based approach and applies K-core decomposition an exact but efficient *O(E+V)* algortihm to identify repetitive regions in metagenomes. We use paired end read information to 
connect unitigs both vertically (same read mapping to multiple unitigs)  and horizontally (paired reads denoting adjacency and unitgs bordering repeats) in order to preserve homology information in the graph. K-core decomposition then
hierarchially decomposes the graph to reveal unitigs grouped together by abundance of repeats (copy number) and visualized by peaks in the KOMB profile graph (Shells vs Number of unitigs in shells). We test KOMB on simulated, synthetic 
and real metagenomic data. More details can be found in the paper. </p>

## Dependencies
**Dependencies:** In order to run KOMB, you will need [Python 3.x](https://www.python.org/download/releases/3.0/),  [ABySS](https://github.com/bcgsc/abyss), [Bowtie2](https://github.com/BenLangmead/bowtie2)
and [C++ -11](http://www.cplusplus.com/).

You also need [igraph-c](https://igraph.org/c/), a popular C library for large scale network analysis.

KOMB has two major dependencies for Core genome analysis:
```
1. Igraph C Library
2. Boost C++ Library
```
**Optional Dependencies:** [Kraken](https://github.com/DerrickWood/kraken) for running in metagenomic mode, [Spades](http://spades.bioinf.spbau.ru/release3.11.1/manual.html) for GFA input/output. Spades has already been included in external.
## Preparing  data for KOMB
**Read filtering:**
We use [kmer_filter](http://catchenlab.life.illinois.edu/stacks/comp/kmer_filter.php) to filter out the reads. This is included in external. The default filtering setting we use is given by the following string
```
external/kmer_filter -1 READ1 -2 READ2 -o output_filtered -D --abundant --k-len 15 --max_k_freq 2
```
This will filter out all reads containing more than 80% of abundant kmers. The kmer size is 15 and abundance threshold is 2 occurences. 
``--abundant`` flag means that those reads that have abundant kmers will be discarded and ``-D`` flag allows us to capture the discarded reads.
We then use the files ``output_filtered/READ*discards.fq`` for the rest of the process.

## Setup
**Setup:** Clone the KOMB repo and run ```make```.
Make will create a build directory with two subdirectories, namely, apps and objects. The komb executable will be in ``/build/apps/``.

## Run KOMB
**Usage:**
``
python3 komb.py [Flags]
``
The description of various flags is given below. Please note that while running with **-m** it is necessary to specifiy both **-k** and **-db** i.e path to kraken and the database to use for taxonomonic classification
```
usage: komb.py [-h] [-m] [-s] [-1 READ1] [-2 READ2] [-c] [-g GENOMESIZE]
              [-l LEVEL] [-k KRAKEN] [-db DATABASE] [-n NUMHITS] [-e KMER]
              [-f] [-u]

KOMB: K-core decomposition on unitig graph

arguments:
  -h, --help            show this help message and exit
  -m, --classification     Reads are classified  first
  -s, --regular         Default KOMB
  -1 READ1, --read1 READ1
                        P.E Read1.fa/P.E Read1.fq
  -2 READ2, --read2 READ2
                        P.E Read2.fa/P.E Read2.fq
  -c, --correction      Read correction required
  -g GENOMESIZE, --genomesize GENOMESIZE
                        Input genome size
  -l LEVEL, --level LEVEL
                        Classification level for kraken (genus or species)
  -k KRAKEN, --kraken KRAKEN
                        path to kraken
  -db DATABASE, --database DATABASE
                        path to kraken database
  -n NUMHITS, --numhits NUMHITS
                        Bowtie2 maximum hits per read
  -e KMER, --kmer KMER  Set kmer size (less than equal to 100)
  -f, --gfa             Build from SPAdes GFA graph
  -u, --unitig-filter   Filter out unitigs below read length

```

## Generating plots

To generate the figures from the paper refer to ```Notebooks/KOMB Figures.ipynb``` and ```Notebooks/HMP.ipynb```. All of the required data is in the ```Data``` folder, so simply using **Run All** option in Jupyter will build all of the figures and save them locally as .png files.

To generate the Shakya graph plot you will need the [graphviz](https://www.graphviz.org/) package to build a graph layout. We also use [ImageMagick](https://imagemagick.org/index.php) to convert the resulting layout from .ps to .png format. We run graphviz sfdp tool as follows:

```
sfdp -Tps Data/shakya_color_conn_after50.dot > shakya_color_conn_after50.ps
```

Note that most of the configuration is done inside the .dot file. In particular we use Prism algorithm from graphviz to resolve node overlaps and we set edges to be drawn transparently to avoid cluttering the picture. Finally, sfdp tool creates a force-directed layout for the graph.

## Contributors
- Code:
    * Advait Balaji
    * Nicolae Sapoval
- Project:
    * Dr. R.A. Leo Elworth
    * Dr. Santiago Segarra
    * Dr. Todd J. Treangen

In case of any issues please open an issue on Gitlab Issues page.
