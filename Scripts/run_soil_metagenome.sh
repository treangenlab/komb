#!/bin/bash

# Takes in one argument the dataset ID.
# Example usage:
#     run_dataset ERR2282016
function run_dataset {
	# Create directory for filtered reads and filter based on 
	# presence of abundant k-mers. We will only use reads that
	# consist of at least 80% of abudant k-mers. We will call
	# a k-mer abundant if it occurs twice or mroe in the entire
	# dataset. 
	# k = 15
	mkdir $1"_filter";
	kmer_filter -1 $1_1.fq -2 $1_2.fq \
				-o $1_filter -D --abundant --k_len 15 --max_k_freq 2;
	cp $1_filter/$1_1.discards.fq .;
	cp $1_filter/$1_2.discards.fq .;
	# Some consmetic trimming of FQ headers for correct parsing 
	# later on in KOMB
	sed 's/|.*//' $1_1.discards.fq > $1_1.fil.fq;
	sed 's/|.*//' $1_2.discards.fq > $1_2.fil.fq;
	# Run KOMB on the filtered reads, note that unitigs will be 
	# filtered out if they are shorter than the read length
	/home/Users/ns58/komb/repo/run.py -s -1 $1_1.fil.fq -2 $1_2.fil.fq \
									  -e 50 > $1.komb.log 2>&1;
	# Copy over the outputs
	cp kcore.txt $1.kcore.txt;
	cp pagerank.txt $1.pagerank.txt;
	cp final.unitigs.fa $1.final.unitigs.fa;
} 

# Run for all experiments specified in the list given as argument
for info in $(cat $1); do
	info_arr=( $info )
	fname=${info_arr[0]}
	folder=${info_arr[1]} 
	cd ${folder};
	run_dataset fname;
	cd ..;
done