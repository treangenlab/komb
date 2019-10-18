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
#include "cmdline.h"
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
#include <cstdint>
using namespace std;

//create aliases
namespace BOOST_C = boost::container;

//typedef
typedef boost::char_separator<char> separator_type;

//Outputs wall-clock time to checkpoint important operations
void printConsole(const string printThis){
    time_t result = time(nullptr);
    cout<<asctime(localtime(&result))<<" "<<printThis<<endl;
}

//Run K-core and PageRank on the graph
void runAlgorithm(igraph_t &graph, string dir){

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

//Store edges and strand information in GFALink class
class GFALink{

	private:
		uint32_t unitig_a,unitig_b;
		string orientation_a,orientation_b;
	public:
		GFALink(){};
		GFALink(uint32_t unitig_a,uint32_t unitig_b, string  orientation_a, string orientation_b);
		string getUnitigEdge();
		uint32_t getFirstUnitig();
		uint32_t getSecondUnitig();
		string getFirstUnitigOrientation();
		string getSecondUnitigOrientation();
};


GFALink::GFALink(uint32_t unitig_a,uint32_t unitig_b, string orientation_a, string orientation_b){
	this->unitig_a = unitig_a;
	this->unitig_b  = unitig_b;
	this->orientation_a  = orientation_a;
	this->orientation_b = orientation_b;	
}

string GFALink::getUnitigEdge(){
	return this->unitig_a+"$"+this->unitig_b;
}

uint32_t GFALink::getFirstUnitig(){
	return this->unitig_a;
}

uint32_t GFALink::getSecondUnitig(){
	return this->unitig_b;
}

string GFALink::getFirstUnitigOrientation(){
	return this->orientation_a;
}

string GFALink::getSecondUnitigOrientation(){
	return this->orientation_b;
}


//Process GFA 
void processGFA(const string dir, int readLength){
	BOOST_C::vector<GFALink> link;
	BOOST_C::map<uint32_t, uint32_t> contig_ids;
	BOOST_C::map<uint32_t,uint32_t>::iterator it1, it2;
	string filename;
	if(dir == ""){
		filename = "output.gfa";
	}
	else{
		filename = dir + "/" + "output.gfa";
	}
	ifstream inf(filename);
	string line;
	uint32_t id = 0;
	while(getline(inf,line)){
		BOOST_C::vector<string> tokens;
		boost::tokenizer<separator_type> tokenizer(line, separator_type("\t"));
		for(auto it = tokenizer.begin();it!=tokenizer.end();++it){
			tokens.push_back(*it);
		}
		if(tokens[0] == "S"){
			uint32_t unitig_length = tokens[2].length();
			if(unitig_length >= readLength){
				contig_ids[stoi(tokens[1])]  = id;
				id++;
			}
		}
		if(tokens[0] == "L"){
			it1 = contig_ids.find(stoi(tokens[1]));
			it2 = contig_ids.find(stoi(tokens[3]));
			if (it1 != contig_ids.end() && it2 != contig_ids.end()){
				GFALink gl(contig_ids[stoi(tokens[1])],contig_ids[stoi(tokens[3])],tokens[2],tokens[4]);
				link.push_back(gl);
			}
		}
	}
	contig_ids.clear();	
	//Make igraph from links
	igraph_t graph;
	igraph_vector_t edges;
	igraph_vector_init(&edges,0);
	for(auto &n : link){
		igraph_vector_push_back(&edges,n.getFirstUnitig());
		igraph_vector_push_back(&edges,n.getSecondUnitig());
	}
	link.clear();
	igraph_create(&graph ,&edges,0,0);
	igraph_vector_destroy(&edges);
	runAlgorithm(graph,dir);
	
}

//Check if unordered set is a subset
bool is_subset_of(const boost::unordered_set<uint32_t> &a, const boost::unordered_set<uint32_t> &b){
	if (a.size() > b.size()){
		return false;
	}
	auto const not_found = b.end();
	for (auto const& element: a){
		if (b.find(element) == not_found){
			return false;
		}
	}
	return true;
}

//Create custom hash for unordered_set<unordered_set<uint_32> >
namespace customset{
	struct cset{

		boost::unordered_set<uint32_t> uset;

	};
	//overload == to compare boost::unordered_set<uint_32>
	bool operator==(const cset &cs1, const cset &cs2){
		return(cs1.uset == cs2.uset);
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
void  readSam(const string filename, BOOST_C::map<string,boost::unordered_set<uint32_t>  > &u_map){
	ifstream inf(filename);
	string line;
	while(getline(inf,line)){
		BOOST_C::vector<string> tokens;
		boost::tokenizer<separator_type> tokenizer(line, separator_type("\t"));
		for(auto it = tokenizer.begin();it!=tokenizer.end();++it){
			tokens.push_back(*it);

		}
		if(tokens.size()>7){
			int pos = tokens[0].find("/");
			u_map[tokens[0].substr(1,pos)].insert(stoi(tokens[2]));
		}	
	}
}

//Create final set unitigs to connect
BOOST_C::vector<BOOST_C::vector<uint32_t> > processDictionary( BOOST_C::map<string,boost::unordered_set<uint32_t> > &u_map1, BOOST_C::map<string,boost::unordered_set<uint32_t> > &u_map2){
	BOOST_C::vector<BOOST_C::vector<uint32_t> > u_vec;
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
		boost::unordered_set<uint32_t> temp_set;
		temp_set.insert(u_map1[n].begin(), u_map1[n].end());
		temp_set.insert(u_map2[n].begin(), u_map2[n].end());
		if(temp_set.size() >= 2){
			customset::cset temp_cset = {temp_set};
			auto it_set = seen_set.find(temp_cset);
			if(it_set==seen_set.end()){
				seen_set.insert(temp_cset);
				BOOST_C::vector<uint32_t> temp(temp_cset.uset.begin(),temp_cset.uset.end());
				u_vec.push_back(temp);
			}
		}
	}
	printConsole("Created unitig set");
	
	//Step 3: Clear memory
	u_map1.clear();
	u_map2.clear();

	//print number of reads to parse	
	printConsole("Number of reads: "+to_string(u_vec.size()));
	
	return u_vec;
}


void createGraph(BOOST_C::vector<BOOST_C::vector<uint32_t> > &u_vec, uint32_t num_unitig, string dir){	
	//set of seen pair - done to avoid simplifying graphs and writing more edges to file than needed
	set<pair<uint32_t,uint32_t> > edge_set;
	FILE *f;
	if(dir == ""){
		f = fopen("edges.txt","w");
	}
	else{
		f = fopen((dir+"/"+"edges.txt").c_str(),"w");
	}
	#pragma omp parallel for
	for(int i = 0; i<u_vec.size();i++){
		BOOST_C::vector<uint32_t> temp = u_vec[i];
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
	runAlgorithm(graph,dir);
}


//main
int main(int argc, char* argv[]){
	
	//parse cmdline
	cmdline::parser parse;
	parse.add<string>("directory",'d',"Directory to store result",false,"");
	parse.add<int>("readlength",'r',"Read Length",false,100);
	parse.add("gfa",'g',"Run GFA build");

	parse.parse_check(argc, argv);
	
	string dir = parse.get<string>("directory");
	int readLength = parse.get<int>("readlength");
	if(parse.exist("gfa")){
		processGFA(dir,readLength);
		exit(0);
	}
	
	string alignment1,alignment2;
	if(dir == ""){
		alignment1 = "alignment1.sam";
		alignment2 = "alignment2.sam";	
	}
	else{
		alignment1 = dir+"/alignment1.sam";
		alignment2 = dir+"/alignment2.sam";
	}


	//create dictionaries to store reads to unitigs mappings
	BOOST_C::map<string,boost::unordered_set<uint32_t> > u_map1;
	BOOST_C::map<string,boost::unordered_set<uint32_t> > u_map2;
	
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
	uint32_t num_unitig = numUnitig(dir);
	
	//prepare read to unitigs dictionary
	BOOST_C::vector<BOOST_C::vector<uint32_t> > u_vec  = processDictionary(u_map1,u_map2);	
	
	
	//build graph
	createGraph(u_vec, num_unitig, dir);	

	//process graph
	processGraph(dir);
		
	return 0;
}
