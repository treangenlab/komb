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

//edges
boost::unordered_set<pair<int,int> > edges;

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
int numUnitig(){
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
BOOST_C::map<int,boost::unordered_set<int> > processDictionary( BOOST_C::map<string,boost::unordered_set<int> > &u_map1, BOOST_C::map<string,boost::unordered_set<int> > &u_map2){
	BOOST_C::map<int,boost::unordered_set<int> > u_map;
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
	int count = 0;
	for(auto &n :  pair_reads){
		boost::unordered_set<int> temp_set;
		temp_set.insert(u_map1[n].begin(), u_map1[n].end());
		temp_set.insert(u_map2[n].begin(), u_map2[n].end());
		if(temp_set.size() >= 2){
			customset::cset temp_cset = {temp_set};
			auto it_set = seen_set.find(temp_cset);
			if(it_set==seen_set.end()){	
				seen_set.insert(temp_cset);
				u_map[count++] = temp_cset.uset;
			}
		}
	}
	printConsole("Created unitig set");
	//Step 3: Clear memory
	u_map1.clear();
	u_map2.clear();

	return u_map;
}


void createGraph(BOOST_C::map<int,boost::unordered_set<int> > &u_map, int num_unitig){	
	//set of seen pairs
	set<pair<int,int> > edge_set;
	FILE *f;
	f = fopen("edges.txt","w");
	printConsole("Number of reads: "+to_string(u_map.size()));
	#pragma omp parallel for
	for(int i = 0; i<u_map.size();i++){
		BOOST_C::vector<int> temp(u_map[i].begin(),u_map[i].end());
		#pragma omp parallel for collapse(2)
		for(int j = 0; j < temp.size();++j){
			for(int k = 0;  k < temp.size();++k){
				if(j-k > 0){
					if(temp[j] > temp[k]){
						auto it = edge_set.find(make_pair(temp[j],temp[k]));
						if (it == edge_set.end()){
							#pragma omp critical
							{
								fprintf(f,"%d\t%d\n",temp[j],temp[k]);
								//outf<<temp[j]<<"\t"<<temp[k]<<endl;
								edge_set.insert(make_pair(temp[j],temp[k]));
							}
						}
					}
					else{
						auto it = edge_set.find(make_pair(temp[k],temp[j]));
						if (it == edge_set.end()){
							#pragma omp critical
							{
								fprintf(f,"%d\t%d\n",temp[k],temp[j]);
								//outf<<temp[k]<<"\t"<<temp[j]<<endl;
								edge_set.insert(make_pair(temp[k],temp[j]));
							}
						}
					}
				}
			}
		}
	}
	fclose(f);
	//clear memory
        u_map.clear();
	edge_set.clear();	
	/*
	printConsole("Begin iterations");
	#pragma omp parallel for
	for(int i=0;i<u_map.size();i++){
	igraph_sparsemat_t spmat,cmat;
        igraph_vector_int_t rid,cid;
        igraph_vector_int_init(&rid,1);
        igraph_vector_int_init(&cid,1);
	igraph_sparsemat_init(&spmat,num_unitig,num_unitig,UINT_MAX);
	igraph_real_t elem;
		BOOST_C::vector<int> temp(u_map[i].begin(),u_map[i].end());
		#pragma omp parallel for collapse(2)
		for(int j = 0; j < temp.size() ; j++){
			for(int k = 0; k < temp.size(); k++){
				if(j-k > 0){
					if(temp[j] > temp[k]){
						igraph_vector_int_set(&rid,0,temp[j]);
				                igraph_vector_int_set(&cid,0,temp[k]);
                				#pragma omp critical
                				{
                        				igraph_sparsemat_compress(&spmat,&cmat);
                        				igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
                        				++elem;
                        				igraph_sparsemat_entry(&spmat,temp[j],temp[k],1);
                				}
					}
					else{
						igraph_vector_int_set(&rid,0,temp[k]);
                				igraph_vector_int_set(&cid,0,temp[j]);
                				#pragma omp critical
                				{
                        				igraph_sparsemat_compress(&spmat,&cmat);
                        				igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
                        				++elem;
                        				igraph_sparsemat_entry(&spmat,temp[k],temp[j],1);
                				}			
					}
				}
			}
		}	
	}
	printConsole("End iterations");
	//clear memory
	u_map.clear();
	igraph_vector_int_destroy(&rid);
	igraph_vector_int_destroy(&cid);

	igraph_t graph;
	igraph_sparsemat_compress(&spmat,&cmat);
	igraph_sparsemat_droptol(&cmat, 2);
	igraph_sparsemat(&graph, &cmat,0);
	igraph_sparsemat_destroy(&spmat);
	igraph_sparsemat_destroy(&cmat);
	igraph_destroy(&graph);
	printConsole("Finished iteration");
	*/

}

void processGraph(){
	//Read graph into memory
	printConsole("Reading Graph");
	igraph_t graph;
        FILE *inpf;
        inpf = fopen("edges.txt", "r"); 
        igraph_read_graph_edgelist(&graph, inpf, 0, 0);
        igraph_vector_t vec;
        igraph_vector_init(&vec,1);
        printConsole("Running K-core");
	igraph_coreness(&graph, &vec, IGRAPH_ALL);
	igraph_vector_destroy(&vec);
	igraph_destroy(&graph);
}


//main
int main(int argc, char* argv[]){
	
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
			readSam("alignment1.sam",u_map1);
			#pragma omp task
			readSam("alignment2.sam",u_map2);
		}

	}
	printConsole("End reading SAM files");
	
	
	//get unitig  length to create the sparse matrix
	int num_unitig = numUnitig();
	
	//prepare read to unitigs dictionary
	BOOST_C::map<int,boost::unordered_set<int> > u_map = processDictionary(u_map1,u_map2);
	
	//build graph
	createGraph(u_map, num_unitig);
	
	//process graph
	processGraph();	

	/*	
	igraph_t graph;
	FILE *inpf;
	inpf = fopen("edges.txt", "r"); 
	igraph_read_graph_edgelist(&graph, inpf, 0, 0);
	igraph_vector_t vec;
	igraph_vector_init(&vec,1);
	igraph_coreness(&graph, &vec, IGRAPH_ALL);
	*/
	
	/*	
	//graph building code   
        igraph_sparsemat_t spmat,cmat;
	igraph_vector_int_t rid,cid;
        igraph_vector_int_init(&rid,1);
        igraph_vector_int_init(&cid,1);
	igraph_sparsemat_init(&spmat,2,2,1);
        //igraph_sparsemat_entry(&spmat,0,0,0);
        //igraph_sparsemat_entry(&spmat,0,1,0);
        //igraph_sparsemat_entry(&spmat,1,0,0);
        //igraph_sparsemat_entry(&spmat,1,1,3);
	//igraph_vector_int_t rid,cid;
	//igraph_vector_int_init(&rid,1);
	//igraph_vector_int_init(&cid,1);
	//igraph_real_t elem;
	//igraph_vector_int_set(&rid,0,1);
	//igraph_vector_int_set(&cid,0,1);
	//igraph_sparsemat_compress(&spmat,&cmat);
	//igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
	//cout<<elem<<endl;	
	//igraph_sparsemat_entry(&spmat,1,1,elem++);		
	//igraph_sparsemat_compress(&spmat,&cmat);
	//igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
	//cout<<elem<<endl;	
	igraph_real_t elem;
	#pragma omp parallel for 
	for(int i=0; i< 100000; ++i){
		//igraph_sparsemat_t cmat;
		igraph_vector_int_set(&rid,0,1);
        	igraph_vector_int_set(&cid,0,1);
		#pragma omp critical
		{
			//igraph_real_t elem;
			igraph_sparsemat_compress(&spmat,&cmat);
			igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
			++elem;
			igraph_sparsemat_entry(&spmat,1,1,1);
		}
	}
	//igraph_sparsemat_t cmat;
	igraph_vector_int_set(&rid,0,1);
	igraph_vector_int_set(&cid,0,1);
	igraph_sparsemat_compress(&spmat,&cmat);
	igraph_sparsemat_index(&cmat,&rid,&cid, nullptr, &elem);
	cout<<elem<<endl;
	*/	


	/*
	igraph_t graph;
	igraph_adjlist_t adj;
	igraph_adjlist_init_empty(&adj,4);
	igraph_vector_int_t *p = igraph_adjlist_get(&adj,0);
	igraph_vector_int_push_back(p,3);
	igraph_vector_int_t *q = igraph_adjlist_get(&adj,3);
	igraph_vector_int_push_back(q,1);
	//igraph_free(p);
	//cout<<igraph_vector_int_size(p)<<endl;
	igraph_adjlist(&graph, &adj, IGRAPH_ALL, 1);	
	igraph_integer_t vcount = igraph_vcount(&graph);
	igraph_integer_t ecount = igraph_ecount(&graph);
	cout<<"Number of vertices: "<<vcount<<endl;
	cout<<"Number of edges: "<<ecount<<endl;
	igraph_adjlist_destroy(&adj);
	igraph_destroy(&graph);
	*/
	
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

	/*
	boost::unordered_set<int> uset;
	for(int i = 0;i<1000;i++){
		uset.insert(i);
	}

	#pragma omp parallel for 
	for (int j = 0; j<1000;++j){
		itertools::combine(uset);
		cout<<"\r"<<j;
	}
	cout<<endl;
	cout<<edges.size()<<endl;
	*/
	/*
	int n = 1000;
	int r = 2;
	BOOST_C::vector<int> res;
	BOOST_C::vector<bool> v(n);
    	fill(v.end() - r, v.end(), true);
    	do {
        	for (int i = 0; i < n; ++i) {
            		if (v[i]) {
                		//cout << (i + 1) << " ";
				res.push_back(i+1);
            		}
        	}
        	//cout << "\n";
    	} while (next_permutation(v.begin(), v.end()));
	*/
	return 0;
}
