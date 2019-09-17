import sys
import os
from collections import defaultdict
import itertools
from itertools import combinations 
import random
import time
import json
from igraph import *

random.seed(7)

def dumpJson(shells,dir=''):
	with open(os.getcwd()+'/'+dir+'/result_shells.json', 'w') as fp:
		json.dump(shells, fp)

def printConsole(message):
	print >> sys.stderr,time.strftime("%c")+': '+ message

def report(dict_k):
	num_nodes = []
	for i in range(k-1):
		if i+1 == k-1:
			num_nodes.append(len(dict_k[k-1]))
		else:
			core1 = dict_k[i+1]
			core2 = dict_k[i+2]
			shell = [node for node in core1 if node not in core2]
			num_nodes.append(len(shell))
	print(num_nodes)
	print(len(num_nodes))

def buildGraph(read2unitigs, dir=''):
	G = snap.TUNGraph.New()
	#nodes_added = set()
	pairs = []
	printConsole('Building Graph')
	for read in read2unitigs:
		pairs = list(combinations(read2unitigs[read],2))
		for pair in pairs:
			a,b = pair
			G.AddEdge2(a,b)
	printConsole('\tGraph has '+str(G.GetNodes())+' nodes and '+ str(G.GetEdges())+ ' edges')
	printConsole('\tSaving edge list in SNAP format')
	snap.SaveEdgeList(G,os.getcwd()+"/"+dir+'/graph.txt')
	return G

def graphSecond(read2unitigs, dir = ''):
	dict_k = defaultdict(list)
	edge_dict = defaultdict(int)
	printConsole('Building Graph')
	for read in read2unitigs:
	#	edges.extend(list(combinations(read2unitigs[read],2)))
		for a,b in combinations(read2unitigs[read],2):
			if a > b:
				edge_dict[(a,b)]+=1
			else:
				edge_dict[(b,a)]+=1
	read2unitigs.clear()
	for k in list(edge_dict):
		if edge_dict[k] < 3:
			del edge_dict[k]
	G = Graph(edge_dict.keys())
	edge_dict.clear()
	#del  edges[:]
	#printConsole('Simplifying Graph')
	#G = G.simplify()
	printConsole('Running K-core')
	shells = G.shell_index()
	for k in xrange(1,max(shells)+1):
		idx = [i for i in xrange(len(shells)) if shells[i]==k]
		dict_k[k] = idx
		subg = G.induced_subgraph(idx)
		subg.write(os.getcwd()+'/'+dir+'/sub_graph_'+str(k)+'.gml','gml')
	dumpJson(dict_k,dir)
	printConsole('Saved shells as Json')
	printConsole('Running PageRank')
	ranks = G.pagerank()
	with open(os.getcwd()+'/'+dir+'/pagerank.txt','w+') as wf:
		for i in xrange(len(ranks)):
			wf.write(str(i)+'\t'+ str(ranks[i])+'\n')
	G.write(os.getcwd()+'/'+dir+'/complete_graph.gml','gml')
	printConsole('Saved complete graph as GML')
	return G


def read_sam(file):
	printConsole('Building Dictionary: '+file.split('/')[-1])
	read2unitigs = defaultdict(set)
	with open(file,'r') as inpf:
		for line in inpf:
			tokens = line.split('\t')
			if len(tokens) > 7:
				read2unitigs[tokens[0][:-2]].add(int(tokens[2]))
	return read2unitigs

def processDictionary(read2unitigs1,read2unitigs2):
	read2unitigs = defaultdict(set)
	seen_node_set = set()
	printConsole('Eliminating single pairs')
	pair_reads = [id for id in read2unitigs1 if id in read2unitigs2]
	printConsole('Preparing final Dictionary')
	printConsole('\tTruncating set of unitigs to 1000')
	for read in pair_reads:
		read2unitigs1[read].update(read2unitigs2[read])
		if len(read2unitigs1[read]) >= 2:
			if tuple(read2unitigs1[read]) in seen_node_set:
				continue
			seen_node_set.add(tuple(read2unitigs1[read]))
			if len(read2unitigs1[read]) > 1000:
				read2unitigs[read] = set(random.sample(read2unitigs1[read],1000))
			else:
				read2unitigs[read] = read2unitigs1[read]
	printConsole('\tFinished discarding reads mapped to same unitig')
	printConsole('\tFinished discarding reads with duplicate unitig sets')
	printConsole('\tFinal number of reads: '+str(len(read2unitigs)))
	printConsole('Clearing Memory')
	read2unitigs1.clear()
	read2unitigs2.clear()
	return read2unitigs

def main():
	cwd = os.getcwd()
	read2unitigs1 = read_sam(cwd+'/alignment1.sam')
	read2unitigs2 = read_sam(cwd+'/alignment2.sam')
	read2unitigs = processDictionary(read2unitigs1,read2unitigs2)
	G = buildGraph(read2unitigs)
	dict_k = runKCore(G)
	report(dict_k)
	runPageRank(G)
	printConsole('Finished')



if __name__== "__main__":
	main()
