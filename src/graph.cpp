//
// Created by Advait Balaji on 6/30/20.
//

#include <cstring>
#include "graph.h"
#include <omp.h>
#include <algorithm>

using namespace gfa;

namespace komb
{

    Kgraph::Kgraph(int threads)
    {
        _threads = threads;
    }

    Kgraph:: Kgraph(int threads, int readlength)
    {
        _threads = threads;
        _readlength = readlength;
    }

    void Kgraph::readSAM(const std::string &samfile, umapset &umap)
    {
        FILE* f = fopen(samfile.c_str(),"r");
        if (f == nullptr) { exit(EXIT_FAILURE);}
        char* line = nullptr;
        size_t len;
        ssize_t read;

        while(read = getline(&line,&len,f) != -1)
        {
            if(line[0] != '@')
            {
                char* token = strtok(line,"\t");
                std::string read(token);
                uint8_t count = 0;
                while(token && (count < 2))
                {
                    token = strtok(NULL, "\t");
                    count++;
                }
                std::string unitig(token);
                umap[read.substr(1,read.find('/'))].insert(std::stoi(unitig));
            }

        }
    }

    vvec Kgraph::getEdgeInfo(umapset &umap1, umapset &umap2)
    {
        vvec edgeinfo;
        usset seen_set;
        std::vector<std::string> pair_read;

        //filtering singleton reads
        for(auto & it : umap1)
        {
            std::string read = it.first;
            if (umap2.find(read) != umap2.end())
            {
                pair_read.emplace_back(read);
            }
        }

        //filter those sets already seen or that are a subset of previously seen sets
        int num_reads = pair_read.size();
        fprintf(stdout,"Num pairs to be processed: %d\n",num_reads);
        //int workload = num_reads/_threads;
        #pragma omp parallel for num_threads(_threads)
        for(uint32_t i = 0; i < num_reads; i++)
        {
            std::set<uint32_t> tempset;
            #pragma omp critical
            {
                tempset.insert(umap1[pair_read[i]].begin(), umap1[pair_read[i]].end());
                tempset.insert(umap2[pair_read[i]].begin(), umap2[pair_read[i]].end());
                if (tempset.size() > 1) 
                {
                    if (seen_set.find(tempset) == seen_set.end()) 
                    {
                        seen_set.insert(tempset);
                        std::vector<uint32_t> tempvec{tempset.begin(), tempset.end()};
                        edgeinfo.emplace_back(tempvec);
                    }
                }
            }
        }
        umapset().swap(umap1);
        umapset().swap(umap2);
        return edgeinfo;
    }



    void Kgraph::generateGraph(vvec& vec, const std::string& dir)
    {
        std::string edgelist_file = dir+"/edgelist.txt";
        uspair seen_edge;
        FILE* ef = fopen(edgelist_file.c_str(),"w+");
        int num_sets_unitigs = vec.size();
        fprintf(stdout,"Number of unitig sets to interconnect: %zu\n",vec.size());
        for(int i = 0; i < num_sets_unitigs; i++)
        {
            std::vector<uint32_t> tempvec = vec[i];

            #pragma omp parallel for collapse(2) num_threads(_threads)
            for(uint32_t j = 0 ; j < tempvec.size(); j++)
            {
                for(uint32_t k = 1; k < tempvec.size(); k++)
                {
                    if(j > k)
                    {
                        //prevent multiedges by only considering a particular ordering of nodes
                        if(tempvec[j] > tempvec[k])
                        {
                            //check if  edge already seen
                            if (seen_edge.find(std::make_pair(tempvec[j],tempvec[k])) == seen_edge.end())
                            {
                                #pragma omp critical
                                {
                                    fprintf(ef,"%d\t%d\n",tempvec[j],tempvec[k]);
                                    seen_edge.insert(std::make_pair(tempvec[j],tempvec[k]));
                                }
                            }
                        }
                        else
                        {
                            //check if edge already seen
                            if (seen_edge.find(std::make_pair(tempvec[k],tempvec[j])) == seen_edge.end())
                            {
                                #pragma omp critical
                                {
                                    fprintf(ef,"%d\t%d\n",tempvec[k],tempvec[j]);
                                    seen_edge.insert(std::make_pair(tempvec[k],tempvec[j]));
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(ef);
        vvec().swap(vec);
        uspair().swap(seen_edge);
    }

    void Kgraph::readEdgeList(const std::string& dir)
    {
        std::string edgelist_file = dir+"/edgelist.txt";
        igraph_t graph;
        FILE* inpf = fopen(edgelist_file.c_str(), "r");
        igraph_read_graph_edgelist(&graph, inpf, 0, 0);
        fclose(inpf);
        fprintf(stdout,"GraphInfo...\n\tNumber of vertices: %d\n",(int)igraph_vcount(&graph));
        fprintf(stdout,"\tNumber of edges: %d\n",(int)igraph_ecount(&graph));

        runCore(graph,dir);
    }

    void Kgraph::runCore(igraph_t &graph, const std::string& dir)
    {
        std::string kcore_file = dir+"/kcore.txt";
        igraph_vector_int_t vec, deg;
        igraph_vector_int_init(&vec,1);
        igraph_vector_int_init(&deg,1);
        igraph_degree(&graph,&deg,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS);
        igraph_coreness(&graph,&vec, IGRAPH_ALL);
        FILE* kcf = fopen(kcore_file.c_str(),"w+");
        int koresize = igraph_vector_int_size(&vec);
        for(int i = 0; i < koresize ; i++)
        {
            fprintf(kcf, "%d\t%d\t%d\n",i,(int)igraph_vector_int_e(&vec,i),(int)igraph_vector_int_e(&deg,i));
        }

        fclose(kcf);
        igraph_vector_int_destroy(&vec);
        igraph_vector_int_destroy(&deg);
        igraph_destroy(&graph);
    }

    void Kgraph::processGFA(const std::string &dir, bool weight)
    {
        std::string gfa_file = dir + "/output.gfa";
        std::vector<Gfa> link;
        std::set<uint32_t> u_id;
        FILE* gfafile = fopen(gfa_file.c_str(),"r");
        if (gfafile == nullptr) { exit(EXIT_FAILURE);}
        char* line = nullptr;
        size_t len;
        ssize_t rec;
        while(rec = getline(&line,&len,gfafile) != -1)
        {
            std::vector<std::string> tokens;
            char* token = strtok(line,"\t");
            std::string id(token); //store the id
            uint8_t count = 0;
            while(token && (count < 4)) // need the first three postions
            {
                std::string token = strtok(NULL, "\t");
                tokens.emplace_back(token);
                count++;
            }
            if(id == "S")
            {
                //is a segment!
                uint32_t ulength = tokens[1].length();
                if(ulength >= _readlength) {
                    u_id.insert(std::stoi(tokens[0]));
                }
            }
            else if(id == "L")
            {
                //is a link!
                auto it1 = u_id.find(std::stoi(tokens[0]));
                auto it2 = u_id.find(std::stoi(tokens[2]));

                if(it1 != u_id.end() && it2 != u_id.end())
                {
                    Gfa gfa_link(std::stoi(tokens[0]),std::stoi(tokens[2]),
                            tokens[1],tokens[3]);
                    link.emplace_back(gfa_link);
                }

            }
            else
            {   //not important!
                continue;
            }
        }
        std::set<uint32_t>().swap(u_id);
        igraph_t graph;
        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges,0);
        for(auto &n : link)
        {
            igraph_vector_int_push_back(&edges,n.getFirstUnitig());
            igraph_vector_int_push_back(&edges,n.getSecondUnitig());
        }
        std::vector<Gfa>().swap(link);
        igraph_create(&graph,&edges,0,0);
        igraph_vector_int_destroy(&edges);
        std::string edgelist_file = dir+"/edgelist.txt";
        FILE* wf = fopen(edgelist_file.c_str(),"w+");
        igraph_write_graph_edgelist(&graph, wf);
        fclose(wf);
        runCore(graph,dir);
        anomalyDetection(dir, weight);
    }
   
    std::map<std::string, std::string> Kgraph::readUnitigsFile(const std::string& dir, bool isBifrost)
    {
        std::map<std::string, std::string> unitigs;
        const std::string utgf = dir + "/unitigs.fasta";
        FILE* fp = fopen(utgf.c_str(),"r");
        if (fp == nullptr) { exit(EXIT_FAILURE);}
        std::string unitig_num; 
        
        char* line = nullptr;
        size_t len;
        ssize_t read;
       
        if (!isBifrost)
        {
            while (read = getline(&line,&len,fp) != -1)
            {
                std::string uline(line); 
                if (uline.substr(0,1) == ">")
                {
                    unitig_num = uline.substr(1,uline.find(" ")-1);
                }
                else
                {
                    unitigs[unitig_num]  = uline.substr(0,uline.length()-1);
                }    
            }    
        }
        else
        {
            while (read = getline(&line,&len,fp) != -1)
            {
                std::string uline(line); 
                if (uline.substr(0,1) == ">")
                {
                    unitig_num = uline.substr(1,uline.length()-2);
                }
                else
                {
                    unitigs[unitig_num]  = uline.substr(0,uline.length()-1);
                }    
            }
        }
        return unitigs;
    }
    
    void Kgraph::combineFile(const std::string& dir, bool isBifrost)
    {    
        std::vector<std::vector<std::string> > kcore_lines;
        std::map<std::string, std::string> unitigs;
        
        const std::string kcf = dir + "/kcore.txt";

        FILE* f = fopen(kcf.c_str(),"r");
        if (f == nullptr) { exit(EXIT_FAILURE);}
        char* line = nullptr;
        size_t len;
        ssize_t read;

        while(read = getline(&line,&len,f) != -1)
        {
            std::vector<std::string> cur_line;
            char* token = strtok(line,"\t");
            while(token)
            {
                std::string cur_token(token);
                cur_line.emplace_back(cur_token);
                token = strtok(NULL, "\t");
            }    
       
            //cur_line[1] = cur_line[1].substr(0,cur_line[1].length()-1);
            cur_line[1].erase(std::remove(cur_line[1].begin(), cur_line[1].end(), '\n'), cur_line[1].end());
            kcore_lines.emplace_back(cur_line);
        }
 
        unitigs = Kgraph::readUnitigsFile(dir,isBifrost);
        
        const std::string cmbf = dir + "/combined.fasta";
    
        FILE* outf = fopen(cmbf.c_str(),"w+");
        for(uint64_t i = 0; i < kcore_lines.size(); i++)
        {   
            std::string header = ">Unitig_"+kcore_lines[i][0]+"|"+ kcore_lines[i][1];
            fprintf(outf, "%s\n", header.c_str());
            auto key = unitigs.find(kcore_lines[i][0]);
            if (key != unitigs.end())
            {
                fprintf(outf, "%s\n", key->second.c_str());
            }
        }
   
        fclose(outf);
    }
    
    void Kgraph::createRER(long long int& vertices, long long int& edges)
    {
        /* test igraph random graph generation */                                
        igraph_t er_rand;                                                  
        igraph_rng_seed(igraph_rng_default(),7);                           
        igraph_erdos_renyi_game(&er_rand, IGRAPH_ERDOS_RENYI_GNM, vertices, edges,0,0); /* bool directed, bool loops */
        igraph_vector_int_t vec;                                               
        igraph_vector_int_init(&vec,1);                                        
        igraph_coreness(&er_rand,&vec, IGRAPH_ALL);                           
        /* delete the above block to remove igraph random graph test */
     }
    
    void Kgraph::anomalyDetection(const std::string& dir, bool weight)
    {
                                
        std::string edgelist_file = dir+"/edgelist.txt";                   
        igraph_t graph;                                                    
        FILE* inpf = fopen(edgelist_file.c_str(), "r");                                                         
        igraph_read_graph_edgelist(&graph, inpf, 0, 0);                                                         
        igraph_lazy_adjlist_t al;                                          
        igraph_lazy_adjlist_init(&graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);                                   
        CombineCoreA::run(al,dir,weight);   
        
    }

    double Kgraph::getMedian(std::vector<double> vec, int start, int end)
    {
       double median = 0;
       int size = end - start - 1;
       if (size % 2 == 0)
       {
           median = (vec[start + size/2 -1] + vec[start + size/2]) / 2;

       }
       else
       {
           median = vec[start + (size - 1) / 2];
       }
       return median;

    }

    void Kgraph::splitAnomalousUnitigs(const std:: string& dir, bool isBifrost)
    {
        std::string anomalyFile = dir+"/top_scoring_anomalous_unitigs.txt";
        std::string backgroundFile = dir+"/low_scoring_anomalous_unitigs.txt";
        std::string coreAf =  dir+"/CoreA_anomaly.txt";
        std::map<std::string, std::string> unitigs = Kgraph::readUnitigsFile(dir, isBifrost);
        std::vector<std::vector<std::string> > coreA_lines;
        FILE* inp_coreAf = fopen(coreAf.c_str(),"r");
     
        if (inp_coreAf == nullptr) { exit(EXIT_FAILURE);}
        char* line = nullptr;
        size_t len;
        ssize_t read;

        while(read = getline(&line,&len,inp_coreAf) != -1)
        {
            std::vector<std::string> cur_line;
            char* token = strtok(line,"\t");
            while(token)
            {
                std::string cur_token(token);
                cur_line.emplace_back(cur_token);
                token = strtok(NULL, "\t");
            }    
       
            //cur_line[1] = cur_line[1].substr(0,cur_line[1].length()-1);
            cur_line[1].erase(std::remove(cur_line[1].begin(), cur_line[1].end(), '\n'), cur_line[1].end());
            coreA_lines.emplace_back(cur_line);
        }
       
        int size = coreA_lines.size();
        std::vector<double> anomaly_scores;
        for(int i = 0; i < size; i++)
        {
            anomaly_scores.emplace_back(std::stod(coreA_lines[i][1]));
        }
        
        std::vector<double> anomaly_scores_copy(anomaly_scores);
        std::sort(anomaly_scores.begin(), anomaly_scores.end());

        double q1 = Kgraph::getMedian(anomaly_scores, 0, size/2 - 1);
        double q3;

        if (size % 2 == 0) 
        {
           q3 = Kgraph::getMedian(anomaly_scores, size/2, size - 1);
        } 
        else 
        {
           q3 = Kgraph::getMedian(anomaly_scores, size/2 + 1, size - 1);
        }
        
        double cutoff = q3 + 1.5 * (q3 - q1); 
        
        FILE* anomalyf = fopen(anomalyFile.c_str(), "w+");
        FILE* backgroundf = fopen(backgroundFile.c_str(), "w+");        

        for (int i = 0; i < size; i++)
        {
            if (anomaly_scores[i] >= cutoff)
            {
                std::string unitig_header = "Unitig_"+std::to_string(i);
                fprintf(anomalyf,"%s\n",unitig_header.c_str());
                auto key = unitigs.find(std::to_string(i));
                if (key != unitigs.end())
                {
                    fprintf(anomalyf,"%s\n", key->second.c_str());
                }
            }
            else
            {
                std::string unitig_header = "Unitig_"+std::to_string(i);
                fprintf(backgroundf,"%s\n",unitig_header.c_str());
                auto key = unitigs.find(std::to_string(i));
                if (key != unitigs.end())
                {
                    fprintf(backgroundf,"%s\n", key->second.c_str());
                }
            }
        }
        
        fclose(anomalyf);
        fclose(backgroundf);
        
    }

} // namespace komb
