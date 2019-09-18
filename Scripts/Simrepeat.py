'''
    Author: Advait Balaji
    Lab: Treangen Lab
    Rice University
'''
import numpy as np
from fastnumbers import fast_int,fast_float
import time
import json
import os,sys
import random

random.seed(7)
#1000, 100, 10 copy_num  for repeats of length .
#min 100 spacer.

#Return the ranks of each element in an integer sample.
def ranks(sample):
    indices = sorted(range(len(sample)), key=lambda i: sample[i])
    return sorted(indices, key=lambda i: indices[i])

#Sample of k elements from range(n), with a minimum distance d.
def sample_with_minimum_distance(n, k, d):
    sample = random.sample(range(n-(k-1)*(d-1)), k)
    return [s + (d-1)*r for s, r in zip(sample, ranks(sample))]


#Create a random sequence
def createSeq(length):
    dict_seq = {1:"A",2:"T",3:"G",4:"C"}
    random_seq = [dict_seq[random.randrange(1, 5)] for _ in range(length)]
    return ''.join(random_seq)

#Mutate certain positions in the sequence
def mutate(seq, positions):

    seq = [char for char in seq]
    for position in positions:
        if seq[position] == "A":
            seq[position] =  random.choice(["T","G","C"])
        elif seq[position] == "G":
            seq[position] = random.choice(["T","A","C"])
        elif seq[position] == "C":
            seq[position] = random.choice(["T","G","A"])
        else:
            seq[position] = random.choice(["A","G","C"])
    return ''.join(seq)

#create positions for min 80% similarity
def createEighty(seq,length, tandem = False):
    eighty = []
    if length == 100:
        copy_num = 999
    elif length == 1000:
        copy_num = 99
    else:
        copy_num = 9
    if tandem:
        copy_num = 19
    for _ in range(copy_num):
        percent =  random.randrange(0,21)
        positions = random.sample(range(0,length), (percent*length)//100)
        eighty.append(mutate(seq,positions))
    return eighty

#create positions for min 90% similarity
def createNinety(seq, length, tandem = False):
    ninety = []
    if length == 100:
        copy_num = 999
    elif length == 1000:
        copy_num = 99
    else:
        copy_num = 9
    if tandem:
        copy_num = 19
    for _ in range(copy_num):
        percent =  random.randrange(0,11)
        positions = random.sample(range(0,length), (percent*length)//100)
        ninety.append(mutate(seq,positions))
    return ninety

#create final dict containing 9 sequences of length 100,1000,10000 with minimum seq identity 80%,90% and 100%
def createFamilyDict(master_seq, lengths):
    dict_family = {}
    for seq,length in zip(master_seq,lengths):
        dict_family[str(length)+"_80"] = [seq]+createEighty(seq,length, False)
        dict_family[str(length)+"_90"] = [seq]+createNinety(seq,length, False)
        if length == 100:
            copy_num = 1000
        elif length == 1000:
            copy_num = 100
        else:
            copy_num = 10
        dict_family[str(length)+"_100"] = [seq for _ in range(copy_num)]
    return dict_family


#create simulgenomes
def createRepeatGenome(dict_family):
    dir = '/Users/advaitbalaji/Downloads/'
    with open(dir+'ecoli.fna','r') as ef:
        lines = ef.readlines()
        ecolihead = lines[0]
        #ecoli_lines = ''.join(ef.readlines()[1:])#.replace('\n','')
        ecoli_lines = ("".join(lines[1:])).replace('\n','')
    #print(ecoli_lines)
    if not os.path.exists(dir+'EcoliForSim_postest'):
        os.mkdir(dir+'EcoliForSim_postest')
    dict_insertionpos = {}
    for key in dict_family.keys():
        #with open(dir+'EcoliForSim/ecoli_'+key+'.fastq','w+') as of:
        copy_num = len(dict_family[key])
        positions = sample_with_minimum_distance(len(ecoli_lines), copy_num, 100)
        with open(dir+'EcoliForSim_postest/ecoli_'+key+'.fastq','w+') as of, open(dir+'EcoliForSim_postest/GroundTruth.txt','a+') as wf:
            ecoli_lines_mod =  ecoli_lines
            assert(len(ecoli_lines_mod) == len(ecoli_lines))
            #print(key,len(positions[:len(dict_family[key])]))
            for seq,pos in zip(dict_family[key],positions[:copy_num]):
                ecoli_lines_mod = ecoli_lines_mod[:pos]+seq+ecoli_lines_mod[pos:]
                wf.write(key+':')
                wf.write(' '+str(pos)+' ')
                wf.write(' '+seq+'\n')
            of.write(ecolihead.split('\n')[0]+'_'+key+'\n')
            of.write(ecoli_lines_mod+'\n')
        print("Written file: ",key)

def createTandem(length, copy_num):
    dict_family = {}
    seq = createSeq(length)
    dict_family[str(length)+"_80"] = [seq]+createEighty(seq,length, True)
    dict_family[str(length)+"_90"] = [seq]+createNinety(seq,length, True)
    dict_family[str(length)+"_100"] = [seq for _ in range(copy_num)]
    with open('ecoli.fna','r') as ef:
        lines = ef.readlines()
        ecolihead = lines[0]
        ecoli_lines = ("".join(lines[1:])).replace('\n','')
    dir = '/Users/advaitbalaji/Downloads/'
    if not os.path.exists(dir+'EcoliTandem'):
        os.mkdir(dir+'EcoliTandem')
    for key in dict_family:
        intra_position = random.sample(range(0,10), copy_num/2)
        inter_position = sample_with_minimum_distance(len(ecoli_lines), copy_num/2 , 30000)
        positions = []
        for i in range(copy_num/2):
            positions.append(inter_position[i])
            positions.append(inter_position[i]+intra_position[i])
        with open(dir+'EcoliTandem/ecoli_'+key+'.fastq','w+') as of, open(dir+'EcoliTandem/GroundTruth.txt','a+') as wf:
            ecoli_lines_mod =  ecoli_lines
            assert(len(ecoli_lines_mod) == len(ecoli_lines))
            count = 0
            for seq,pos in zip(dict_family[key],positions):
                ecoli_lines_mod = ecoli_lines_mod[:pos]+seq+ecoli_lines_mod[pos:]
                wf.write(key+':')
                wf.write(' '+str(pos)+' ')
                wf.write(' '+seq+'\n')
            of.write(ecolihead.split('\n')[0]+'_'+key+'\n')
            of.write(ecoli_lines_mod+'\n')
        print("Written file: ",key)

    ''' 

        ecoli_lines = ("".join(lines[1:])).replace('\n','')
    if not os.path.exists(dir+'EcoliForSim_postest'):
        os.mkdir(dir+'EcoliForSim_postest')
    dict_insertionpos = {}
    '''

#main
def main():
    #lengths = [100,1000,10000]
    #master_seq = [createSeq(fast_int(lengths[0])), createSeq(fast_int(lengths[1])), createSeq(fast_int(lengths[2]))]
    #dict_family = createFamilyDict(master_seq,lengths)
    #createRepeatGenome(dict_family)
    createTandem(10000,20)
        
if __name__== "__main__":
    main()


    #print(master_seq)
