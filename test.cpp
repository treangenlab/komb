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

#define SIZE 900000
using namespace std;

//create aliases
namespace BOOST_C = boost::container;

//typedef
typedef boost::char_separator<char> separator_type;

//create combinations
namespace itertools{
	class combination{
		public:
			//TODO
			void combine(int offset, int k){
				if(k==0){
			//TODO
				
				}	
			

			}

	};

}

//Create custom hash for unordered_set<unordered_set<int> >
namespace customset{
	struct cset{

		boost::unordered_set<int> uset;

	};

	//overload == to compare boost::unordered_set<int>
	bool operator==(const cset &cs1, const cset &cs2){
		bool flag_size = false;
		bool flag_elem = true;
		if(cs1.uset.size() == cs2.uset.size()){
			flag_size = true;
			for(auto &n : cs1.uset){
				auto it = cs2.uset.find(n);
				if (it == cs2.uset.end()){
					flag_elem == false;
					break;
				}
			}
			return flag_size && flag_elem;
		}
		else{
			return false;
		}
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
int unitigLength(){
	FILE *f = NULL;
	int num_lines;
	string command = "wc -l final.unitigs.fa";
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
BOOST_C::map<string,boost::unordered_set<int> > processDictionary( BOOST_C::map<string,boost::unordered_set<int> > &u_map1, BOOST_C::map<string,boost::unordered_set<int> > &u_map2){
	BOOST_C::map<string,boost::unordered_set<int> > u_map;
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
	//Step 2: Create set of unitigs to connect for each pair in Set 1 
	for(auto &n : pair_reads){
		boost::unordered_set<int> temp_set;
		temp_set.insert(u_map1[n].begin(), u_map1[n].end());
		temp_set.insert(u_map[n].begin(),u_map2[n].end());
		customset::cset temp_cset = {temp_set};
		auto it_set = seen_set.find(temp_cset);
		if(it_set!=seen_set.end()){
			continue;
		}
		else{
			seen_set.insert(temp_cset);
			u_map[n] = temp_cset.uset;
		}
	}
	printConsole("Created unitig set");
	//Step 3: Clear memory
	u_map1.clear();
	u_map2.clear();

	return u_map;
}


igraph_t createGraph(BOOST_C::map<string,boost::unordered_set<int> > &u_map){



	//TODO



}


//Main
int main(int argc, char* argv[]){
	//create dictionaries to store reads to unitigs mappings
	BOOST_C::map<string,boost::unordered_set<int> > u_map1;
	BOOST_C::map<string,boost::unordered_set<int> > u_map2;
	
	//test readSam function
	printConsole("Begin reading SAM files");
	#pragma omp parallel
	{
		#pragma omp single
		{
			#pragma omp task
			readSam("alignment1.sam",u_map1);
			#pragma omp task
			readSam("alignment2.sam",u_map2);
		}

	}
	printConsole("End reading SAM files");
	
	
	//get unitig  length to create the sparse matrix
	int unitig_length = unitigLength();
	
	//prepare read to unitigs dictionary
	BOOST_C::map<string,boost::unordered_set<int> > u_map = processDictionary(u_map1,u_map2);

	/*
	//use the igraph C library for operations
	igraph_t graph;
	igraph_sparsemat_t spmat;
	igraph_sparsemat_t cmat;
	igraph_sparsemat_init(&spmat,3,3,2);
	igraph_sparsemat_entry(&spmat,0,0,2);
	igraph_sparsemat_entry(&spmat,0,1,0);
	igraph_sparsemat_entry(&spmat,1,0,0);
	igraph_sparsemat_entry(&spmat,0,2,0);
	igraph_sparsemat_entry(&spmat,1,1,3);
	igraph_sparsemat_entry(&spmat,1,2,0);
	igraph_sparsemat_entry(&spmat,2,0,3);
	igraph_sparsemat_entry(&spmat,2,1,3);
	igraph_sparsemat_entry(&spmat,2,2,3);
	igraph_sparsemat_compress(&spmat,&cmat);
	igraph_sparsemat_destroy(&spmat); 
	igraph_sparsemat_droptol(&cmat, 2);
	igraph_sparsemat(&graph,&cmat,0);
	igraph_integer_t vcount = igraph_vcount(&graph);
	igraph_integer_t ecount = igraph_ecount(&graph);
	cout<<"Number of vertices: "<<vcount<<endl;
	cout<<"Number of edges: "<<ecount<<endl;
	igraph_sparsemat_destroy(&cmat);
	igraph_destroy(&graph);
	cout<<"No Error :)"<<endl;
	*/
	return 0;
}
