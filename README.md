# KOMB - K(c)ore genome analyzer

KOMB uses K-Core and PageRank to analyze core genomes and repeats in single and metagenomes

**Dependencies:** 2. In order to run KOMB, you will need [Python 2.7.x](https://www.python.org/downloads/), [Kraken](https://github.com/DerrickWood/kraken),  [Bcalm2](https://github.com/GATB/bcalm)
and [Lighter](https://github.com/mourisl/Lighter).

You also need _igraph-c_, a popular C library for large scale network analysis as described [here](https://igraph.org/c/).

KOMB has two major dependencies for Core genome analysis:
```
1. Igraph C Library
2. Boost C++ Library
```
**Usage:**
``
python run.py [Flags]
``
The description of various flags is given below. Please note that while running with **-m** it is necessary to specifiy both **-k** and **-db** i.e path to kraken and the database to use for taxonomonic classification
```
usage: run.py [-h] [-m] [-s] [-1 READ1] [-2 READ2] [-c] [-g GENOMESIZE]
              [-l LEVEL] [-k KRAKEN] [-db DATABASE]

Kore Genome Analyzer: Graph based analysis

optional arguments:
  -h, --help            show this help message and exit
  -m, --metagenome      Reads are metagenomes
  -s, --single          Reads are single genomes/related genomes
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

```
You can use the C++ igraph script as well if required by tweaking some functions. Please run
```g++ -O3 test.cpp -Iinclude/igraph/ -Llib/ -ligraph -I/include/boost/ -fopenmp -o test```
