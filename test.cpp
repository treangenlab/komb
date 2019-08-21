extern "C"
{
#include <igraph.h>
}
#include <boost/container/vector.hpp>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>


#define SIZE 900000
using namespace std;

void printConsole(string printThis){
    time_t result = time(nullptr);
    cout<<asctime(localtime(&result))<<" "<<printThis<<endl;
}


int main(){
	boost::container::vector<int>v1;
	v1.push_back(7);
	igraph_t graph;
	igraph_sparsemat_t mat;
	igraph_vector_t cores;
	igraph_sparsemat_init(&mat, SIZE, SIZE, 4294967295);
	long i;
	printConsole("Start populating adjmat");
	for (i=0;i<SIZE-200000;i++){
		igraph_sparsemat_entry(&mat,i,i+1,1);
	}
	printConsole("Built adjmat");
	igraph_sparsemat(&graph,&mat,0);
	printConsole("Built graph");
	igraph_sparsemat_destroy(&mat);
	printConsole("Start Kcore");
	igraph_vector_init(&cores,SIZE);
	igraph_coreness(&graph,&cores, IGRAPH_ALL);
	printConsole("Finished Kcore");
	for(i=0;i<SIZE;i++){
		cout<<igraph_vector_e(&cores,i)<<endl;
	}
	igraph_vector_destroy(&cores);
	igraph_destroy(&graph);
	return 0;

}
