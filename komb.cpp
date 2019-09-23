/*
 *Author: Advait Balaji
 *Date: August 2019
 *
 *Use igraph-C library to run K-Core and PageRank analysis
 */
extern "C"
{
#include <igraph.h>
#include <limits.h>
}
#include <boost/container/map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/container/vector.hpp>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <utility>
#include <set>
#define SIZE 900000
using namespace std;

//create aliases
namespace BOOST_C = boost::container;

//typedef
typedef boost::char_separator<char> separator_type;


namespace  itertools{
	void combine(boost::unordered_set<int> &uset){
		BOOST_C::vector<int> temp(uset.begin(),uset.end());
		omp_set_num_threads(omp_get_max_threads());
		#pragma omp parallel for collapse(2) 
		for(int i = 0; i< temp.size(); ++i){
			for(int j = 0; j< temp.size();++j){
				if(i-j > 0){
					#pragma omp critical
					{
						edges.insert(make_pair(i,j));
					}
				}
			}
		}
	}
}	


//Create custom hash for unordered_set<unordered_set<int> >
namespace customset{
	struct cset{

		boost::unordered_set<int> uset;

	};
	//overload == to compare boost::unordered_set<int>
	bool operator==(const cset &cs1, const cset &cs2){
		if(cs1.uset.size() <= cs2.uset.size()){
			for(auto &n : cs1.uset){
				auto it = cs2.uset.find(n);
				if (it == cs2.uset.end()){
					return false;
				}
			}
		}
		else{
			return false;
		}
		return true;
	}
	//create hash value for boost::unordered_set<int> 
	size_t hash_value(const cset &cs){
		boost::hash<int> hasher;
		int seed;
		int sum = 0;
		for(auto &n : cs.uset){
			sum += n;
		}
		seed = sum+cs.uset.size();
		return hasher(seed);
	}
}

//Outputs wall-clock time to checkpoint important operations
void printConsole(const string printThis){
    time_t result = time(nullptr);
    cout<<asctime(localtime(&result))<<" "<<printThis<<endl;
}

//execute wc -l to count number of unitigs
int numUnitig(string dir){
	FILE *f = NULL;
	int num_lines;
	string command;
	if(dir == ""){
		command = "wc -l final.unitigs.fa";
	}
	else{
		command = "wc -l "+dir+"/final.unitigs.fa";
	}
	f = popen(command.c_str(), "r");
	int err = fscanf(f, "%d", &num_lines);
	pclose(f);
	cout<<"Number of unitigs: "<<to_string(num_lines/2)<<endl;
	return num_lines/2;
}

//read SAM file
void  readSam(const string filename, BOOST_C::map<string,boost::unordered_set<int>  > &u_map){
	ifstream inf(filename);
	string line;
	while(getline(inf,line)){
		BOOST_C::vector<string> tokens;
		boost::tokenizer<separator_type> tokenizer(line, separator_type("\t"));
		for(auto it = tokenizer.begin();it!=tokenizer.end();++it){
			tokens.push_back(*it);

		}
		if(tokens.size()>7){
			u_map[tokens[0].substr(1,tokens[0].length()-2)].insert(stoi(tokens[2]));
		}	
	}
}

//Create final set unitigs to connect
BOOST_C::vector<BOOST_C::vector<int> > processDictionary( BOOST_C::map<string,boost::unordered_set<int> > &u_map1, BOOST_C::map<string,boost::unordered_set<int> > &u_map2){
	BOOST_C::vector<BOOST_C::vector<int> > u_vec;
	boost::unordered_set<customset::cset> seen_set;
	BOOST_C::vector<string> pair_reads;
	printConsole("Filter single pairs");
	//Step 1: Filter reads with pairs
	for(auto it = u_map1.begin();it!=u_map1.end();++it){
		string read = it->first;
		auto it_map = u_map2.find(read);
		if(it_map!=u_map2.end()){
			pair_reads.push_back(read);	
		}
	}
	printConsole("Finished filtering single pairs");
	//Step 2: Create set of unitigs to connect for each pair in set 1 
	for(auto &n :  pair_reads){
		boost::unordered_set<int> temp_set;
		temp_set.insert(u_map1[n].begin(), u_map1[n].end());
		temp_set.insert(u_map2[n].begin(), u_map2[n].end());
		if(temp_set.size() >= 2){
			customset::cset temp_cset = {temp_set};
			auto it_set = seen_set.find(temp_cset);
			if(it_set==seen_set.end()){	
				seen_set.insert(temp_cset);
				BOOST_C::vector<int> temp(temp_cset.uset.begin(),temp_cset.uset.end());
				u_vec.push_back(temp);
			}
		}
	}
	printConsole("Created unitig set");
	//Step 3: Clear memory
	u_map1.clear();
	u_map2.clear();

	return u_vec;
}


void createGraph(BOOST_C::vector<BOOST_C::vector<int> > &u_vec, int num_unitig, string dir){	
	//set of seen pair - done to avoid simplifying graphs and writing more edges to file than needed
	set<pair<int,int> > edge_set;
	FILE *f;
	if(dir == ""){
		f = fopen("edges.txt","w");
	}
	else{
		f = fopen((dir+"/"+"edges.txt").c_str(),"w");
	}
	printConsole("Number of reads: "+to_string(u_vec.size()));
	#pragma omp parallel for
	for(int i = 0; i<u_vec.size();i++){
		BOOST_C::vector<int> temp = u_vec[i];
		#pragma omp parallel for collapse(2)
		for(int j = 0; j < temp.size();++j){
			for(int k = 0;  k < temp.size();++k){
				//Only need to account for unique combinations
				if(j-k > 0){
					//prevent multiedges by only considering a particular ordering of nodes
					if(temp[j] > temp[k]){
						//check if  edge already seen
						auto it = edge_set.find(make_pair(temp[j],temp[k]));
						if (it == edge_set.end()){
							#pragma omp critical
							{
								fprintf(f,"%d\t%d\n",temp[j],temp[k]);
								edge_set.insert(make_pair(temp[j],temp[k]));
							}
						}
					}
					else{
						//check if edge already seen 
						auto it = edge_set.find(make_pair(temp[k],temp[j]));
						if (it == edge_set.end()){
							#pragma omp critical
							{
								fprintf(f,"%d\t%d\n",temp[k],temp[j]);
								edge_set.insert(make_pair(temp[k],temp[j]));
							}
						}
					}
				}
			}
		}
	}
	//close file
	fclose(f);
	
	//clear memory
        u_vec.clear();
	edge_set.clear();	
}

void processGraph(string dir){
	//Read graph into memory
	printConsole("Reading Graph");
	igraph_t graph;
        FILE *inpf;
        if(dir == ""){
		inpf = fopen("edges.txt", "r");
	}
	else{
		inpf = fopen((dir+"/"+"edges.txt").c_str(),"r");
	}	
        igraph_read_graph_edgelist(&graph, inpf, 0, 0);
        fclose(inpf);

	//K-Core implementation
	igraph_vector_t vec;
        igraph_vector_init(&vec,1);
        printConsole("Running K-core");
	igraph_coreness(&graph, &vec, IGRAPH_ALL);
	FILE *outf;
        if(dir == ""){
		outf = fopen("kcore.txt", "w");
	}
	else{
		outf = fopen((dir+"/"+"kcore.txt").c_str(),"w");
	}	
	for(int i = 0; i < igraph_vector_size(&vec); ++i){
		fprintf(outf, "%d\t%d\n", i,int(igraph_vector_e(&vec, i)));
	}
	fclose(outf);
	igraph_vector_destroy(&vec);
	
	//PageRank implementation from PRPACK and on undirected graph
	igraph_vector_t pgvec;
	igraph_vector_init(&pgvec,1);
        printConsole("Running PageRank");
	igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &pgvec, 0, igraph_vss_all(), 0, 0.31, 0, 0);
	FILE *pgf;
        if(dir == ""){
		pgf = fopen("pagerank.txt", "w");
	}
	else{
		pgf = fopen((dir+"/"+"pagerank.txt").c_str(),"w");
	}	
	for(int i = 0; i < igraph_vector_size(&pgvec); ++i){
		fprintf(pgf, "%d\t%.3f\n", i, float(igraph_vector_e(&pgvec, i)));
	}
	fclose(pgf);
	igraph_vector_destroy(&pgvec);
	
	//Destroy Graph
	igraph_destroy(&graph);
}


//main
int main(int argc, char* argv[]){
	
	//Dir is not empty if runnning in metagenome mode
	string dir;
	string alignment1, alignment2;
	if (argc == 2){
		dir = argv[1];
		alignment1  = dir+"/"+"alignment1.sam";
		alignment2  = dir+"/"+"alignment2.sam";
	}
	else{
		dir = "";
		alignment1 = "alignment1.sam";
		alignment2 = "alignment2.sam";
	}


	
	//create dictionaries to store reads to unitigs mappings
	BOOST_C::map<string,boost::unordered_set<int> > u_map1;
	BOOST_C::map<string,boost::unordered_set<int> > u_map2;
	
	//test readSam function
	printConsole("Begin reading SAM files");
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel
	{
		#pragma omp single
		{
			#pragma omp task
			readSam(alignment1,u_map1);
			#pragma omp task
			readSam(alignment2,u_map2);
		}

	}
	printConsole("End reading SAM files");
	
	
	//get unitig  length to create the sparse matrix
	int num_unitig = numUnitig(dir);
	
	//prepare read to unitigs dictionary
	BOOST_C::vector<BOOST_C::vector<int> > u_vec  = processDictionary(u_map1,u_map2);
	
	//build graph
	createGraph(u_vec, num_unitig, dir);
	
	//process graph
	processGraph(dir);	
	return 0;
}
