#!/bin/bash

# Takes in one argument, name of the directory

function run_filtered_reads_unitigs {
	mkdir $1"_filter";
	kmer_filter -1 $1*1*.fastq -2 $1*2*.fastq \
		    -o $1_filter -D --abundant --k_len 21 --max_k_freq 2;
	cp $1_filter/$1*1*discards.fq $1_1.fil.fq;
	cp $1_filter/$1*2*discards.fq $1_2.fil.fq;

	# Run KOMB on the filtered reads, note that unitigs will be
	# filtered out if they are shorter than the read length
	/home/Users/ns58/komb/repo/run.py -s -1 $1_1.fil.fq -2 $1_2.fil.fq \
					  -e 50 -u > $1_fil_uni.komb.log 2>&1;
	# Copy over the outputs
	cp kcore.txt $1_fil_uni.kcore.txt;
	cp pagerank.txt $1_fil_uni.pagerank.txt;
	cp final.unitigs.fa $1_fil_uni.final.unitigs.fa;
}

function run_filtered_reads_nouni {
        # Run KOMB on the filtered reads, note that unitigs WON'T be
        # filtered out
        /home/Users/ns58/komb/repo/run.py -s -1 $1_1.fil.fq -2 $1_2.fil.fq \
                                          -e 50 > $1_fil.komb.log 2>&1;
        # Copy over the outputs
        cp kcore.txt $1_fil.kcore.txt;
        cp pagerank.txt $1_fil.pagerank.txt;
        cp final.unitigs.fa $1_fil.final.unitigs.fa;
}

function run_unfil_reads_unitigs {
        cp $1_filter/$1*1*.fastq $1_1.fq;
        cp $1_filter/$1*2*.fastq $1_2.fq;

        # Run KOMB on the raw reads, note that unitigs will be
        # filtered out if they are shorter than the read length
        /home/Users/ns58/komb/repo/run.py -s -1 $1_1.fq -2 $1_2.fq \
                                          -e 50 -u > $1.komb.log 2>&1;
        # Copy over the outputs
        cp kcore.txt $1.kcore.txt;
        cp pagerank.txt $1.pagerank.txt;
        cp final.unitigs.fa $1.final.unitigs.fa;
}

function run_unfil_reads_nouni {
        # Run KOMB on the raw reads, note that unitigs WON'T be
        # filtered out 
        /home/Users/ns58/komb/repo/run.py -s -1 $1_1.fq -2 $1_2.fq \
                                          -e 50 > $1_nouni.komb.log 2>&1;
        # Copy over the outputs
        cp kcore.txt $1_nouni.kcore.txt;
        cp pagerank.txt $1_nouni.pagerank.txt;
        cp final.unitigs.fa $1_nouni.final.unitigs.fa;
}

cd $1;
echo "Starting run on "$1;
run_filtered_reads_unitigs $1;
run_filtered_reads_nouni $1;
run_unfil_reads_unitigs $1;
run_unfil_reads_nouni $1;
echo "Done."
cd ..;
