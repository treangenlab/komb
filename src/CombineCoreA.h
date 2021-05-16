#ifndef KOMB_COMBINECOREA_H
#define KOMB_COMBINECOREA_H

#include "CoreA.h"
#include "HashIndexedMinHeap.h"
#include <igraph.h>
#include <set>
#include <vector>
#include <climits>
#include <algorithm>

class CombineCoreA
{

	public:

		static void run(igraph_lazy_adjlist_t& al, const std::string& outdir, bool weighted)
		{
			double* nodeSuspiciousness = nullptr;
			CoreA core_a;
			const std::string kcf = outdir+"/kcore.txt";
			std::pair<std::vector<int>, std::vector<int> > coreness_and_degree = core_a.readKOMBOutput(kcf);
			std::vector<int> coreness = coreness_and_degree.first;
			std::vector<int> degree = coreness_and_degree.second;
			double dense_ratio = *std::max_element(coreness.begin(), coreness.end()) / 2;
			fprintf(stdout,"Dense Ratio: %f\n", dense_ratio);
			int n = degree.size();

			if(weighted)
			{
				nodeSuspiciousness = core_a.getAnomalyScore(degree, coreness);
				double max_dmp = *std::max_element(nodeSuspiciousness, nodeSuspiciousness+n);
				fprintf(stdout,"Max CoreA score: %f\n", max_dmp);
				double theoretical_weight = dense_ratio  * (1 / max_dmp);
				//fprintf(stdout,"Theoretical weight: %f\n", theoretical_weight);
				std::string anomaly_output = outdir+"/CoreA_anomaly.txt";
				FILE* fp = fopen(anomaly_output.c_str(), "w+");
				for(int i = 0; i < n; i++)
				{
					fprintf(fp,"%d\t%f\n",i,nodeSuspiciousness[i]);
				}
				fclose(fp);

            	    //for(int i=0; i < degree.size(); i++) 
            	    //{
                    //	nodeSuspiciousness[i] = theoretical_weight * nodeSuspiciousness[i];
            	    //}
			}
            /***
			CombineCoreA cca;
			std::pair<std::vector<int>, std::vector<int> > result = cca.runMerge(al, n, nodeSuspiciousness);
			std::set<int> anomalies;
            
            		for(auto node : result.first) 
            		{
                		anomalies.insert(node);
            		}

            		for(auto node : result.second) 
            		{
                		anomalies.insert(node);
            		}

            		std::string weighted_output = outdir +"/weighted_anomaly_output.txt";
        		FILE* fp = fopen(weighted_output.c_str(),"w+");
        		fprintf(fp, "vertex_index\n");
        	
        		for(auto node : anomalies) 
        		{
        			fprintf(fp,"%d\n",node);
        		}

        		fclose(fp);
             ***/

        		free(nodeSuspiciousness);

		}

		std::pair<std::vector<int>, std::vector<int> > runMerge(igraph_lazy_adjlist_t& al, int& numNodes, double* suspiciousness)
		{
			double* rowDegree = nullptr; //need to free
			double* colDegree = nullptr; //need to free

			double suspiciousSum = 0;

			if(suspiciousness == nullptr)
			{
				rowDegree = (double*)calloc(numNodes,sizeof(double));
				colDegree = (double*)calloc(numNodes,sizeof(double));
			}
			else
			{

				rowDegree = (double*)calloc(numNodes,sizeof(double));
				colDegree = (double*)calloc(numNodes,sizeof(double));
				memcpy(rowDegree,suspiciousness,sizeof(double)*numNodes);
				memcpy(colDegree,suspiciousness,sizeof(double)*numNodes);
				for(int i = 0; i < numNodes; i++)
				{
					suspiciousSum += 2*suspiciousness[i];
				}
			}

			long edgeNum = 0;
			int size;
			for(int src = 0; src < numNodes; src++)
			{
				igraph_vector_t* vec = igraph_lazy_adjlist_get(&al,src);
				size = igraph_vector_size(vec);
				if(size)
				{
					for(int i = 0; i < size; i++)
					{
						rowDegree[src] += 1;
						colDegree[(int)igraph_vector_e(vec,i)] += 1;
						edgeNum += 1;
					}
				}
			}
			
			suspiciousSum += edgeNum;

			HashIndexedMinHeap rowHeap(numNodes);
			for(int i=0; i<numNodes; i++) 
			{
            			rowHeap.insert(i, rowDegree[i]);
        		}

        		HashIndexedMinHeap colHeap(numNodes);
        		for(int j=0; j<numNodes; j++) 
        		{
            			colHeap.insert(j, colDegree[j]);
        		}

        		int* modes = (int*)malloc(sizeof(int)*(2*numNodes)); //need to free
        		int* order = (int*)malloc(sizeof(int)*(2*numNodes)); //need to free

        		bool* removed[2]; // need to free!
        		for(int i = 0; i < 2; i++)
        		{
        			removed[i] = (bool*)malloc(sizeof(bool)*numNodes);
        		}


        		double maxDensity = 0; // density of the densest block so far
        		int maxDensityNodesNum = 0; // number of nodes belongs to the densest block
        		int numOfNodesBelong = 2*numNodes; // number of nodes belongs to the current block
        		while(numOfNodesBelong >= 1) 
        		{

            			std::pair<int,double> rowPair = rowHeap.peek();
            			std::pair<int,double> colPair = colHeap.peek();

            			std::pair<int,double> pair = std::make_pair(INT_MIN,INT_MIN);
            			int modeToRemove = 0;

            			if(((rowPair.first != INT_MIN) && (rowPair.second != INT_MIN)) && 
				   (((colPair.first == INT_MIN) && (colPair.second == INT_MIN))|| rowPair.second < colPair.second)) 
            			{
                			pair = rowHeap.poll();
                			modeToRemove = 0;
            			}
            			else 
            			{
                			pair = colHeap.poll();
                			modeToRemove = 1;
            			}

            			suspiciousSum -= pair.second;
            			int node = pair.first;
            			order[--numOfNodesBelong] = node;
            			modes[numOfNodesBelong] = modeToRemove;

            			double density = suspiciousSum/numOfNodesBelong;
            			if (numOfNodesBelong >= 1 && density > maxDensity) 
            			{
                			maxDensity = density;
                			maxDensityNodesNum = numOfNodesBelong;
            			}

            			removed[modeToRemove][node] = true;
            
            			int size;
           			igraph_vector_t* vec = igraph_lazy_adjlist_get(&al,node);
    				size = igraph_vector_size(vec);
            			if(modeToRemove == 0) 
            			{ //remove from rows
                			for (int i = 0; i < size; i++) 
              				{
                    				int dst = igraph_vector_e(vec,i);
                    				if (!removed[1][dst]) 
                    				{
                        				colHeap.refreshPriority(dst, colHeap.getPriority(dst) - 1);
                    				}
                			}
            			}
            			else 
            			{ //remove from columns
            				for (int i = 0; i < size; i++) 
                			{
                    				int dst = igraph_vector_e(vec,i);
                    				if (!removed[0][dst]) 
                    				{
                        				rowHeap.refreshPriority(dst, rowHeap.getPriority(dst) - 1);
                    				}
                			}
            			}
        		}

        		int numRowsBelongs = 0;
        		int numColsBelongs = 0;
        		for(int i=0; i<maxDensityNodesNum; i++) 
        		{
            			if(modes[i] == 0) 
            			{
                			numRowsBelongs++;
            			}
            			else 
            			{
                			numColsBelongs++;
            			}
        		}

        		std::vector<int> rows(numRowsBelongs, 0); 
        		std::vector<int> cols(numRowsBelongs, 0); 
        	
        		int rowIndex = 0;
        		int colIndex = 0;
        	
        		for(int i=0; i<maxDensityNodesNum; i++) 
        		{
            			if(modes[i] == 0) 
            			{
            	    			rows[rowIndex++] = order[i];
            			}
            			else 
            			{
                			cols[colIndex++] = order[i];
            			}
        		}

        		for(int i = 0; i < 2; i++)
        		{
        			free(removed[i]);
        		}

        		free(rowDegree);
        		free(colDegree);
        		free(modes);
        		free(order);

        		return std::make_pair(rows,cols);
    		}

};
#endif
