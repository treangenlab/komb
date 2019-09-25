import os
import argparse
import sys
import time
import subprocess
from subprocess import Popen, PIPE
import kraken
from kraken import *
from collections import defaultdict

def cmd_exists(cmd):
	return subprocess.call("type " + cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def getReadLength(file):
	count = 0
	with open(file,"r") as fp:
		for i, line in enumerate(fp):
			if i == 1:
				count = len(line)
				break
	return count

def callLighter(genomesize,read1,read2,concat=True):
	print(time.strftime("%c")+': Starting read correction',file=sys.stderr)
	if not cmd_exists('lighter'):
		print(time.strftime("%c")+': Lighter does not exist in PATH (or not installed), exiting process',file=sys.stderr)
		sys.exit(1)
	print(time.strftime("%c")+': Using size for genome as '+ str(genomesize)+' base pairs',file=sys.stderr)
	print(time.strftime("%c")+': Running error correction using Lighter',file=sys.stderr)
	try:
		p = subprocess.check_output('lighter -K ' + " 32 "+ str(genomesize) + " -t 100 "+ "-r "+ read1 + " -maxcor 20 -noQual", shell=True)
		print(time.strftime("%c")+': Finished Read 1',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error correcting read 1',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('lighter -K ' + " 32 " +str(genomesize) + " -t 100 "+ "-r "+ read2 + " -maxcor 20 -noQual", shell=True)
		print(sys.stderr,time.strftime("%c")+': Finished Read 2',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error correcting read 2',file=sys.stderr)
		sys.exit(1)
	if concat:
		try:
			p = subprocess.check_output('cat  *.cor.* > final.fa ', shell=True)
			print(time.strftime("%c")+': Concatenated reads for Bcalm',file=sys.stderr)
		except subprocess.CalledProcessError as err:
			print(time.strftime("%c")+': Error concatenating reads',file=sys.stderr)
			sys.exit(1)

def callBcalmMeta(folder):
	try:
		p = subprocess.check_output('cat '+folder+'/*.fa > '+ folder+'/final.fa ', shell=True)
		print(time.strftime("%c")+': Concatenated reads for Bcalm',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error concatenating reads',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('bcalm -in '+folder+'/final.fa -kmer-size 33 -abundance-min 5' , shell=True)
		print(time.strftime("%c")+': Unitigs created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error creating unitigs',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('rm -rf *.unitigs.fa.*', shell=True)
		print(time.strftime("%c")+': Temp bcalm files deleted',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error deleting temp bcalm files',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('rm -rf *.h5', shell=True)
		print(time.strftime("%c")+': Temp .h5 file deleted',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error deleting temp .h5 file',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('mv final.unitigs.fa '+folder+'/', shell=True)
		print(time.strftime("%c")+': Moving final unitigs to folder',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error moving unitigs',file=sys.stderr)
		sys.exit(1)


def callBowtie2Meta(read1size,read2size,folder,numhits):
	try:
		p = subprocess.check_output('bowtie2-build '+folder+'/final.unitigs.fa idx', shell=True)
		print(time.strftime("%c")+': Bowtie2 index created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error creating index',file=sys.stderr)
		sys.exit(1)	
	try:
		p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+'  -3 '+str(read1size-33)+' -U '+ folder+'/reads1.fa'+' -p 100 --no-unal > '+folder+'/alignment1.sam', shell=True)
		print(time.strftime("%c")+': Alignment1.sam created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error in bowtie2-align -> Alignment1',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+' -5 '+str(read2size-33)+' -U '+folder+ '/reads2.fa' +' -p 100 --no-unal > '+ folder+'/alignment2.sam', shell=True)
		print(time.strftime("%c")+': Alignment2.sam created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error in bowtie2-align -> Alignment2',file=sys.stderr)
		sys.exit(1)	
	try:
		p = subprocess.check_output('rm -rf idx* ', shell=True)
		print(time.strftime("%c")+': Removed index files',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error removing index files',file=sys.stderr)
		sys.exit(1)

def callBcalm():
	try:
		p = subprocess.check_output('bcalm -in final.fa -kmer-size 33 -abundance-min 5', shell=True)
		print(time.strftime("%c")+': Unitigs created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error creating unitigs',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('rm -rf *.unitigs.fa.*', shell=True)
		print(time.strftime("%c")+': Temp bcalm files deleted',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error deleting temp bcalm files',file=sys.stderr)
		sys.exit(1)
	try:
		p = subprocess.check_output('rm -rf *.h5', shell=True)
		print(time.strftime("%c")+': Temp .h5 file deleted',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error deleting temp .h5 file',file=sys.stderr)
		sys.exit(1)

def callBowtie2(read1size,readfile1,read2size,readfile2,ext1,ext2,mode,numhits):
	try:
		p = subprocess.check_output('bowtie2-build final.unitigs.fa idx', shell=True)
		print(time.strftime("%c")+': Bowtie2 index created',file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error creating index',file=sys.stderr)
		sys.exit(1)
	if mode == 'SC':
		if ext1 == '.fq':
			try:
				p = subprocess.check_output('bowtie2 -x idx --sensitive -k '+str(numhits)+' -3 '+str(read1size-33)+' -U '+ readfile1+'.cor'+ext1+' -p 100 --no-unal > alignment1.sam', shell=True)
				print(time.strftime("%c")+': Alignment1.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment1',file=sys.stderr)
				sys.exit(1)
			try:
				p = subprocess.check_output('bowtie2 -x idx --sensitive -k '+str(numhits)+' -5 '+str(read2size-33)+' -U '+ readfile2+'.cor'+ext2+' -p 100 --no-unal > alignment2.sam', shell=True)
				print(time.strftime("%c")+': Alignment2.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment2',file=sys.stderr)
				sys.exit(1)
		else:
			try:
				p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+' -3 '+str(read1size-33)+' -U '+ readfile1+'.cor'+ext1+' -p 100 --no-unal > alignment1.sam', shell=True)
				print(time.strftime("%c")+': Alignment1.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment1',file=sys.stderr)
				sys.exit(1)
			try:
				p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+' -5 '+str(read2size-33)+' -U '+ readfile2+'.cor'+ext2+' -p 100 --no-unal > alignment2.sam', shell=True)
				print(time.strftime("%c")+': Alignment2.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment2',file=sys.stderr)
				sys.exit(1)
	elif mode == 'SNC':
		if ext1 == '.fq':
			try:
				p = subprocess.check_output('bowtie2 -x idx --sensitive -k '+str(numhits)+' -3 '+str(read1size-33)+' -U '+ readfile1+ext1+' -p 100 --no-unal > alignment1.sam', shell=True)
				print(time.strftime("%c")+': Alignment1.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment1',file=sys.stderr)
				sys.exit(1)
			try:
				p = subprocess.check_output('bowtie2 -x idx --sensitive -k '+str(numhits)+' -5 '+str(read2size-33)+' -U '+ readfile2+ext2+' -p 100 --no-unal > alignment2.sam', shell=True)
				print(time.strftime("%c")+': Alignment2.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment2',file=sys.stderr)
				sys.exit(1)
		else:
			try:
				p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+' -3 '+str(read1size-33)+' -U '+ readfile1+ext1+' -p 100 --no-unal > alignment1.sam', shell=True)
				print(time.strftime("%c")+': Alignment1.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment1',file=sys.stderr)
				sys.exit(1)
			try:
				p = subprocess.check_output('bowtie2 -f -x idx --sensitive -k '+str(numhits)+' -5 '+str(read2size-33)+' -U '+ readfile2+ext2+' -p 100 --no-unal > alignment2.sam', shell=True)
				print(time.strftime("%c")+': Alignment2.sam created',file=sys.stderr)
			except subprocess.CalledProcessError as err:
				print(time.strftime("%c")+': Error in bowtie2-align -> Alignment2',file=sys.stderr)
				sys.exit(1)
	else:
		print(time.strftime("%c")+': Mode unrecognized, exiting process',file=sys.stderr)
		sys.exit(1)
		try:
			p = subprocess.check_output('rm -rf idx* ', shell=True)
			print(time.strftime("%c")+': Removed index files',file=sys.stderr)
		except subprocess.CalledProcessError as err:
			print(time.strftime("%c")+': Error removing index files',file=sys.stderr)
			sys.exit(1)

def callMetagenomePipeline(correction,genomesize,read1,read2,level,kraken,database,numhits):
	mode = 'M'
	read1size = 0
	read2size = 0
	readfile1 = ''
	readfile2 = ''
	ext1 =''
	ext2 = ''
	read1_cor = ''
	read2_cor = ''
	num_reads = 0

	read1size = int(getReadLength(read1))
	readfile1 = read1.split('/')[-1][:-3]
	ext1 = read1.split('/')[-1][-3:]

	read2size = int(getReadLength(read2))
	readfile2 = read2.split('/')[-1][:-3]
	ext2 = read2.split('/')[-1][-3:]	
	try:
		assert(ext1 == ext2)
	except AssertionError as ar:
		print(time.strftime("%c")+': Both read files need to have same extention',file=sys.stderr)
		sys.exit(1)

	try:
		p = subprocess.check_output('wc -l '+read1,shell = True)
		if ext1 == '.fq':
			num_reads = int(p.split(' ')[0])/4
		else:
			num_reads = int(p.split(' ')[0])/2
		print(time.strftime("%c")+': Number of reads: '+ str(num_reads),file=sys.stderr)
	except subprocess.CalledProcessError as err:
		print(sys.stderr,time.strftime("%c")+': Error counting reads',file=sys.stderr)
	if correction:
		callLighter(genomesize,read1,read2,False)
		read1_cor = readfile1.split(ext1)[0]+'.cor'+ext1
		read2_cor = readfile2.split(ext1)[0]+'.cor'+ext2
	else:
		read1_cor = read1
		read2_cor = read2


	k_obj = Kraken(read1_cor,read2_cor,level,kraken,database,num_reads)
	k_obj.runKraken()
	print(time.strftime("%c")+': Started reading taxonomy files',file=sys.stderr)
	k_obj.readNodes(database+'/taxonomy/nodes.dmp')
	k_obj.readNames(database+'/taxonomy/names.dmp')
	temp_dict = k_obj.readKraken()
	filtered_dict = k_obj.filterReadLevel(temp_dict)
	dir_keys = k_obj.separateReads(filtered_dict)

	for key in dir_keys:
		callBcalmMeta('SORT_'+str(key))
		callBowtie2Meta(read1size,read2size,'SORT_'+str(key),numhits)
		#read2unitigs1 = kga.read_sam('SORT_'+str(key)+'/alignment1.sam')
		#read2unitigs2 = kga.read_sam('SORT_'+str(key)+'/alignment2.sam')
		#read2unitigs = kga.processDictionary(read2unitigs1,read2unitigs2)
		#G = kga.graphSecond(read2unitigs,'/SORT_'+str(key))
		try:
			p = subprocess.check_output('./komb '+'SORT_'+str(key))
			print(p.strip('\n'))
		except subprocess.CalledProcessErroe as err:
			print(time.strftime("%c")+': Error running KOMB',file=sys.stderr)
		print(sys.stderr,time.strftime("%c")+': Finished',file=sys.stderr)

def callSinglegenomePipeline(correction,genomesize,read1,read2,numhits):
	mode = ''
	read1size = 0
	read2size = 0
	readfile1 = ''
	readfile2 = ''
	ext1 =''
	ext2 = ''
	num_reads = 0

	read1size = int(getReadLength(read1))
	readfile1 = read1.split('/')[-1][:-3]
	ext1 = read1.split('/')[-1][-3:]

	read2size = int(getReadLength(read2))
	readfile2 = read2.split('/')[-1][:-3]
	ext2 = read2.split('/')[-1][-3:]

	try:
		assert(ext1 == ext2)
	except AssertionError as ar:
		print(time.strftime("%c")+': Both read files need to have same extention',file=sys.stderr)
		sys.exit(1)

	if correction:
		callLighter(genomesize,read1,read2,True)
		mode = 'SC'
	else:
		try:
			p = subprocess.check_output('cat  '+ read1 + ' ' + read2+ ' > final.fa ', shell=True)
			print(time.strftime("%c")+': Concatenated reads for Bcalm',file=sys.stderr)
			mode = 'SNC'
		except subprocess.CalledProcessError as err:
			print(time.strftime("%c")+': Error concatenating reads',file=sys.stderr)
			sys.exit(1)

	callBcalm()
	callBowtie2(read1size,readfile1,read2size,readfile2,ext1,ext2,mode,numhits)
	#cwd = os.getcwd()
	#read2unitigs1 = kga.read_sam(cwd+'/alignment1.sam')
	#read2unitigs2 = kga.read_sam(cwd+'/alignment2.sam')
	#read2unitigs = kga.processDictionary(read2unitigs1,read2unitigs2)
	#G = kga.graphSecond(read2unitigs)
	try:
		p = subprocess.check_output('./komb')
		print(p.strip('\n'))
	except subprocess.CalledProcessError as err:
		print(time.strftime("%c")+': Error running KOMB',file=sys.stderr)
	print(time.strftime("%c")+': Finished',file=sys.stderr)


def main():
	cwd=os.path.dirname(os.path.abspath(__file__))
	parser = argparse.ArgumentParser(description="Kore Genome Analyzer: Graph based analysis")
	parser.add_argument("-m","--metagenome", help="Reads are metagenomes", action = 'store_true')
	parser.add_argument("-s","--single",help="Reads are single genomes", action = 'store_true')
	parser.add_argument("-1",'--read1',type = str, help="P.E Read1.fa/P.E Read1.fq")
	parser.add_argument("-2","--read2",type = str, help="P.E Read2.fa/P.E Read2.fq",)
	parser.add_argument("-c","--correction",help="Read correction required",action = 'store_true')
	parser.add_argument("-g","--genomesize", type = int, help="Input genome size", default=200000000)
	parser.add_argument("-l","--level", type=str, help="Classification level for kraken (genus or species)", default='genus')
	parser.add_argument("-k","--kraken", type=str, help="path to kraken")
	parser.add_argument("-db","--database", type=str, help="path to kraken database")
	parser.add_argument("-n","--numhits",type=int, help="Bowtie2 maximum hits per read", default = 1000)
	args = parser.parse_args()

	if args.metagenome and args.single:
		print(time.strftime("%c")+': Both -m and -s are set, exiting process',file=sys.stderr)
		sys.exit(1)
	if not args.metagenome and not args.single:
		print(time.strftime("%c")+': Please set -m OR -s flag depending on your data, exiting process',file=sys.stderr)
		sys.exit(1)
	
	if not args.read1:
		print(time.strftime("%c")+': Read 1 not provided, exiting process',file=sys.stderr)
		sys.exit(1)
	
	if not args.read2:
		print(time.strftime("%c")+': Read 2 not provided, exiting process',file=sys.stderr)
		sys.exit(1)

	if not cmd_exists('bcalm'):
		print(time.strftime("%c")+': Bcalm does not exist in PATH (or not installed), exiting process',file=sys.stderr)
		sys.exit(1)

	if not cmd_exists('bowtie2'):
		print(time.strftime("%c")+': Bowtie2 does not exist in PATH (or not installed), exiting process',file=sys.stderr)
		sys.exit(1)

	if args.metagenome:
		if not args.kraken:
			print(time.strftime("%c")+': Kraken path not given, exiting process',file=sys.stderr)
			sys.exit(1)
		if not args.database:
			print(time.strftime("%c")+': Kraken database path not given, exiting process',file=sys.stderr)
			sys.exit(1)
		print(time.strftime("%c")+': Starting Kore Genome Analysis on metagenome',file=sys.stderr)
		if args.level.lower() =='genus':
			print(time.strftime("%c")+': Kraken output will be grouped by genus',file=sys.stderr)
		elif args.level.lower() == 'species':
			print(time.strftime("%c")+': Kraken output will be grouped by species',file=sys.stderr)
		else:
			print(time.strftime("%c")+': Unidentified level (-l) option: [DEFAULT] Kraken will be grouped by genus',file=sys.stderr)
		callMetagenomePipeline(args.correction,args.genomesize,args.read1,args.read2,args.level.lower(),args.kraken,args.database,args.numhits)

	elif args.single:
		print(time.strftime("%c")+': Starting Kore Genome Analysis on single genome',file=sys.stderr)
		callSinglegenomePipeline(args.correction,args.genomesize,args.read1,args.read2,args.numhits)
	else:
		print(time.strftime("%c")+': Exiting process',file=sys.stderr)
		sys.exit(1)


if __name__ == '__main__':
	main()
