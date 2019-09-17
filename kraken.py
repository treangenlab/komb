import os
from collections import defaultdict
import pickle
import sys
import time
import subprocess
from subprocess import Popen, PIPE


class Kraken:

	def __init__(self,read1,read2,level,kraken,database,num_reads):
		self.read1 = read1
		self.read2 = read2
		self.level = level
		self.kraken = kraken
		self.database = database
		self.num_reads = num_reads

	def readNodes(self,file):
		print >> sys.stderr, time.strftime("%c")+': Reading taxonomy nodes file'
		species_genus = defaultdict(str)
		with open(file,'r') as rf:
			for line in rf:
				tokens = line.split('\t')
				if tokens[4] == 'species':
					species_genus[tokens[0]] = tokens[2]

		with open('tax_sg_file.pkl', 'wb') as output:
			pickle.dump(species_genus, output)

	def readNames(self,file):
		print >> sys.stderr, time.strftime("%c")+': Reading taxonomy names file'
		names = defaultdict(str)
		with open(file,'r') as rf:
			for line in rf:
				tokens = line.split('\t')
				if tokens[6] == 'scientific name':
					names[tokens[0]] = tokens[2]
		with open('tax_names_file.pkl', 'wb') as out:
			pickle.dump(names, out)

	def runKraken(self):
		try:
			p = subprocess.check_output(self.kraken+'/./kraken --preload --db '+ self.database +' --threads 80 --paired '+ self.read1 + ' '+ self.read2+ ' > sequences.kraken ', shell=True)
			print >> sys.stderr, time.strftime("%c")+': Classification complete'
		except subprocess.CalledProcessError as err:
			print >> sys.stderr,time.strftime("%c")+': Error running kraken'
			sys.exit(1)
		

	def readKraken(self):
		print >> sys.stderr, time.strftime("%c")+': Reading Kraken output'
		reads_dict = defaultdict(list)
		s_time = time.time()
		with open(os.getcwd()+'/sequences.kraken','r') as inpf:
			for line in inpf:
				tokens = line.split('\t')
				if tokens[0] == 'C':
					reads_dict[tokens[2]].append(tokens[1])
		return reads_dict

	def filterReadLevel(self,tax_dict):
		print >> sys.stderr, time.strftime("%c")+': Filtering reads according to selected taxonomic level'
		species_genus = defaultdict(str)
		names = defaultdict(str)
		with open('tax_sg_file.pkl', 'rb') as handle:
			species_genus = pickle.load(handle)
		if self.level == 'species':
			for key in list(tax_dict):
				if not key in species_genus:
					del tax_dict[key]
		else:
			for key in list(tax_dict):
				if key in species_genus:
					tax_dict[species_genus[key]].extend(tax_dict[key])
					del tax_dict[key]

		species_genus.clear()
		print >> sys.stderr, time.strftime("%c")+': Removing taxa with less than '+ str(round(0.005*self.num_reads))+' reads'
		for key in list(tax_dict):
			if len(tax_dict[key]) < 0.005*self.num_reads:
				del tax_dict[key]

		if len(tax_dict) == 0:
			print >> sys.stderr, time.strftime("%c")+': Kraken classification failed'
			sys.exit(1)

		print >> sys.stderr, time.strftime("%c")+': Creating file with names of classified taxa in the sample'
		with open('tax_names_file.pkl', 'rb') as handle:
			names = pickle.load(handle)
		with open('Tax_names.txt','w+') as wf:
			for k in tax_dict:
				wf.write(str(k)+'\t'+names[k]+'\n')

		print >> sys.stderr, time.strftime("%c")+': Creating directory for each taxon'
		for k in tax_dict:
			if not os.path.exists('SORT_'+k):
				os.mkdir('SORT_'+k)

		return tax_dict

	def separateReads(self,filtered_dict):

		if '.fq' in self.read1:
			read1_dict = defaultdict(str)
			read2_dict = defaultdict(str)
			print >> sys.stderr, time.strftime("%c")+': Sorting reads into different directories'
			with open(self.read1,'r') as handle1:
				count = 0
				key  = ''
				flag = False
				for line in handle1:
					if count%4 == 0:
						key  = line.split(' ')[0][1:].strip('\n')
						flag = True
						count+=1
					else:
						if flag:
							read1_dict[key]  = line
							count+=1
							flag = False
						else:
							count+=1
			print >> sys.stderr, time.strftime("%c")+': Prepared  dictionary 1'
			writ_count = 0
			for k in filtered_dict:
				with open('SORT_'+str(k)+'/reads1.fa','w+') as wf1:
					for r_id in filtered_dict[k]:
						if r_id in read1_dict:
							wf1.writelines(['>'+r_id+'\n',read1_dict[r_id]])
				writ_count+=1
				print >> sys.stderr, time.strftime("%c")+': Finished processing '+str(writ_count)+'/'+str(len(filtered_dict))

			read1_dict.clear()
			print >> sys.stderr, time.strftime("%c")+': Cleared memory'

			with open(self.read2,'r') as handle2:
				count = 0
				key  = ''
				flag = False
				for line in handle2:
					if count%4 == 0:
						key  = line.split(' ')[0][1:].strip('\n')
						flag = True
						count+=1
					else:
						if flag:
							read1_dict[key]  = line
							count+=1
							flag = False
						else:
							count+=1
			print >> sys.stderr, time.strftime("%c")+': Prepared  dictionary 2'
			writ_count = 0
			for k in filtered_dict:
				with open('SORT_'+str(k)+'/reads2.fa','w+') as wf1:
					for r_id in filtered_dict[k]:
						if r_id in read1_dict:
							wf1.writelines(['>'+r_id+'\n',read1_dict[r_id]])
				writ_count+=1
				print >> sys.stderr, time.strftime("%c")+': Finished processing '+str(writ_count)+'/'+str(len(filtered_dict))
			print >> sys.stderr, time.strftime("%c")+': Finished read sorting'
		
		else:
			read1_dict = defaultdict(str)
			read2_dict = defaultdict(str)
			print >> sys.stderr, time.strftime("%c")+': Sorting reads into different directories'
			with open(self.read1,'r') as handle1:
				count = 0
				key  = ''
				for line in handle1:
					if count%2 == 0:
						key  = line.split(' ')[0][1:].strip('\n')
						count+=1
					else:
						read1_dict[key] = line
						count+=1
			print >> sys.stderr, time.strftime("%c")+': Prepared  dictionary 1'
			writ_count = 0
			for k in filtered_dict:
				with open('SORT_'+str(k)+'/reads1.fa','w+') as wf1:
					for r_id in filtered_dict[k]:
						if r_id in read1_dict:
							wf1.writelines(['>'+r_id+'\n',read1_dict[r_id]])
				writ_count+=1
				print >> sys.stderr, time.strftime("%c")+': Finished processing '+str(writ_count)+'/'+str(len(filtered_dict))

			read1_dict.clear()
			print >> sys.stderr, time.strftime("%c")+': Cleared memory'

			with open(self.read2,'r') as handle2:
				count = 0
				key  = ''
				for line in handle2:
					if count%2 == 0:
						key  = line.split(' ')[0][1:].strip('\n')
						count+=1
					else:
						read1_dict[key] = line
						count+=1
			print >> sys.stderr, time.strftime("%c")+': Prepared  dictionary 2'
			writ_count = 0
			for k in filtered_dict:
				with open('SORT_'+str(k)+'/reads2.fa','w+') as wf1:
					for r_id in filtered_dict[k]:
						if r_id in read1_dict:
							wf1.writelines(['>'+r_id+'\n',read1_dict[r_id]])
				writ_count+=1
				print >> sys.stderr, time.strftime("%c")+': Finished processing '+str(writ_count)+'/'+str(len(filtered_dict))
			print >> sys.stderr, time.strftime("%c")+': Finished read sorting'
		return filtered_dict.keys()


	def checkLength(self, filtered_dict):
		lengths = []
		for key in filtered_dict:
			lengths.append(len(filtered_dict[key]))
		print max(lengths)
		print min(lengths)
		print sum(lengths)/len(lengths)


def main():
	k_obj= Kraken('SRR606249_1.cor.fq','SRR606249_2.cor.fq','genus')
	print >> sys.stderr, time.strftime("%c")+': Started reading taxonomy files'
	k_obj.readNodes('nodes.dmp')
	k_obj.readNames('names.dmp')
	temp_dict= k_obj.readKraken()
	filtered_dict = k_obj.filterReadLevel(temp_dict)
	k_obj.separateReads(filtered_dict)



if __name__ == '__main__':
	main()









