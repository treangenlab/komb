#ifndef KOMB_COREA_H
#define KOMB_COREA_H


#include <vector>
#include <cmath>  
#include <string>
#include <cstring>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cstdlib>


class CoreA
{

	public:
		
		std::pair<std::vector<int>, std::vector<int> > readKOMBOutput(const std::string& kcf)
 		{                                                                      
                                                                       
        		std::vector<int> coreness;                                         
        		std::vector<int> degree;                                           
                                                                         
      			FILE* fp = fopen(kcf.c_str(),"r");                                 
      			if (fp == nullptr) { exit(EXIT_FAILURE);}                          
      			std::string unitig_num;                                            
      			char* line = nullptr;                                              
      			size_t len;                                                        
      			ssize_t read;                                                      
      			while (read = getline(&line,&len,fp) != -1)                        
      			{                                                                  
          			char* token = strtok(line,"\t"); //ignore                      
          			int i = 0;                                                     
          			while (token != NULL)                                          
          			{                                                                                              
              				if(i == 2)                                                 
              				{                                                          
                  				std::string cur_token(token); // convert to string to easily remove "\n"
                 				cur_token.erase(std::remove(cur_token.begin(), cur_token.end(), '\n'), cur_token.end());
                  				degree.emplace_back(stoi(cur_token));                  
                  				i++;                                                   
             				}                                                          
              				if(i == 1)                                                 
              				{                                                          
                  				coreness.emplace_back(atoi(token)); // can directly take kcore
                 				i++;                                                   
              				}
					if(i == 0)
					{
						i++;	
					}
					token = strtok(NULL,"\t"); 
          			}                                                              
      			}                                                                  
                                                                         
      			fclose(fp);                                                        
      			//fprintf(stdout,"Size of coreness:%lu\n",coreness.size());          
      			//fprintf(stdout,"Size of degree:%lu\n",degree.size());              
      			return std::make_pair(coreness, degree);                           
 		}
	
	
		void run(const std::string& kcf, const std::string& output)
		{

			std::pair<std::vector<int>, std::vector<int> > coreness_and_degree = this->readKOMBOutput(kcf);
			std::vector<int> coreness = coreness_and_degree.first;
			std::vector<int> degree = coreness_and_degree.second;

			int n = degree.size();

			double* anomaly_standard = this->getAnomalyScore(degree, coreness);
			double* anomaly_avg = (double*)malloc(sizeof(double)*n); //create a copy
			double* anomaly = (double*)malloc(sizeof(double)*n); //create a copy
			
			memcpy(anomaly_avg, anomaly_standard, sizeof(double)*n); //copy anomaly_standard to different anomaly_avg
			memcpy(anomaly, anomaly_avg, sizeof(double)*n); //copy anomaly_avg to different anomaly
        		double* anomalyRankStandard = this->standardRank(anomaly_standard, n); //here anomaly_standard is freed!

        		int* orderedIndices = (int*)malloc(sizeof(int)*n);
        		for(int i = 0; i < n; i++) 
        		{
            			orderedIndices[(int)anomalyRankStandard[i]-1] = i;
        		}
        	
        		free(anomalyRankStandard); //freed anomalyRankStandard

        		double* anomalyRankAvg = this->fractionalRank(anomaly_avg,n); //here anomaly_avg is freed!

        		FILE* fp = fopen("anomaly_nodes.txt","w+");
        		fprintf(fp,"rank\tvertex_index\tanomaly_score\tcoreness\tdegree\n");
        	
        		for(int i = n - 1; i >= 0; i--) 
        		{
            			int index = orderedIndices[i];
            			std::string line = std::to_string((n-anomalyRankAvg[index]+1)) + "\t" + std::to_string(index) + "\t" 
            			+ std::to_string(anomaly[index]) + "\t" + std::to_string(coreness[index]) + "\t" + std::to_string(degree[index]);
            			fprintf(fp,"%s\n",line.c_str());
        		}

        		fclose(fp); //close file


        		free(orderedIndices); //freed ordered indices
        		free(anomalyRankAvg); //freed anomalyRankAvg
        		free(anomaly); //freed anomaly


		}



		double* getAnomalyScore(std::vector<int>& degree, std::vector<int>& coreness)
		{
			int n = degree.size();
			double* corenessWithDegree = (double*)malloc(sizeof(double)*n);
			double* degree_ptr_arr = (double*)malloc(sizeof(double)*n);

			for(int i = 0; i < n; i++)
			{
				degree_ptr_arr[i] = degree[i];
			}

			for(int i = 0; i < n; i++)
			{
				corenessWithDegree[i] = coreness[i] * n + degree[i];
			}

			double* corenessRank = this->fractionalRank(corenessWithDegree,n);
			double* degreeRank = this->fractionalRank(degree_ptr_arr,n);
			double* anomaly = (double*)malloc(sizeof(double)*n);
			
			for(int i = 0; i < n; i++)
			{
				anomaly[i] =  std::abs(std::log(degreeRank[i]) - std::log(corenessRank[i]));
			}

			//free(corenessWithDegree);
			//free(degree_ptr_arr);
			free(corenessRank);
			free(degreeRank);

			return anomaly;
		}

		double* fractionalRank(double* scores, int& n) 
		{
			double* fractionalRank = (double*)malloc(sizeof(double)*n);
	    		std::vector<int> list;

	    		for (int i = 0; i < n; i++) 
	    		{
	        		list.push_back(scores[i]);
	    		}
	    	
	    		std::sort(list.begin(), list.end(), std::greater<int>{});

	    		list.erase(std::unique(list.begin(), list.end()), list.end());
	 
	    		int rank = 0;

	  		for (auto value : list) 
	    		{
	        		double avg = 0.0;
	        		int cnt = 0;
	 
	        		//for (auto& e : scores) 
	        		for(int i = 0; i < n; i++)
	        		{
	            			if (scores[i] == value) 
	            			{
	                			rank++;
	                			cnt++;
	                			avg += rank;
	            			}
	        		}
	        
	        		avg /= cnt;
	 
	        		for(int i = 0; i < n; i++)
	        		{
	            			if (scores[i] == value) 
	            			{
	            				fractionalRank[i] = avg;
	            			}
	        		}
	    		}

	    		free(scores); // free scores here!
			return fractionalRank;
		}

		double* standardRank(double* scores, int n)
		{

			double* standardRank = (double*)malloc(sizeof(double)*n);
    			std::vector<int> list;
    			for(int i = 0; i < n; i++) 
    			{
        			list.push_back(scores[i]);
    			}

    			std::sort(list.begin(), list.end(), std::greater<int>{});
    			list.erase(std::unique(list.begin(), list.end()), list.end());
 
    			int rank = 1;
    			for (auto value : list) 
    			{
        			int temp = rank;
        			for (int i = 0; i < n; i++) 
        			{
            				if (scores[i] == value) 
            				{
            					standardRank[i] = temp;
                				//std::cout << temp << " " << value << " " << e.first.c_str() << std::endl;
                				rank++;
            				}
        			}
    			}
    			free(scores); //free scores
    			return standardRank;
		}



	private:

		double* transform(int* values, int& n)
		{
			double* result = (double*)malloc(sizeof(double)*n);
			for(int i = 0; i < n; i++)
			{
				result[i] = values[i];
			}
			return result; // result will be freed in fractionalRank
		}
};
#endif
