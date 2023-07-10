//
// Created by Advait Balaji on 06/30/2020.
// Modified by Nicolae Sapoval on 06/25/2023.
// Modified by Marko Tanevski on 07/07/2023
//

#include <cstring>
#include "graph.h"
#include <omp.h>
#include <algorithm>
#include <unordered_set>

using namespace gfa;

namespace komb
{
    Kgraph::Kgraph(uint32_t threads)
    {
        _threads = threads;
    }

    Kgraph::Kgraph(uint32_t threads, uint64_t readlength)
    {
        _threads = threads;
        _readlength = readlength;
    }

    void Kgraph::fileNotFoundError(const std::string& path) {
        std::cerr << "File " << path << " could not be opened. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    void Kgraph::readSAM(const std::string &samfile, umapset &umap, bool fulgor)
    {
        FILE* f = fopen(samfile.c_str(), "r");
        if (f == nullptr) { fileNotFoundError(samfile); }
        char* line = nullptr;
        size_t len;
        ssize_t bytes_read;

        while(bytes_read = getline(&line, &len, f) != -1)
        {
            if (!fulgor) {
                if(line[0] != '@')
                {
                    char* token = strtok(line, "\t");
                    std::string read(token);
                    uint8_t count = 0;
                    while(token && (count < 2))
                    {
                        token = strtok(NULL, "\t");
                        count++;
                    }
                    std::string unitig(token);
                    if (unitig.compare("*") != 0) {  // '*' indicates an unmapped read in SAM
                        umap[read.substr(1, read.find('/'))].insert(unitig);
                    }
                }
            } else {
                char* token = strtok(line, "\t");
                std::string read(token);
                token = strtok(NULL, "\t");
                uint32_t count = std::stoi(token);
                for(uint32_t i = 0; i < count; i++) {
                    token = strtok(NULL, "\t");
                    std::string unitig_id(token);
                    unitig_id.erase(std::remove(unitig_id.begin(), unitig_id.end(), '\n'), unitig_id.cend());
                    umap[read.substr(0, read.find('/'))].insert(unitig_id);
                }
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
        fprintf(stdout, "Num pairs to be processed: %d\n", num_reads);
        //int workload = num_reads/_threads;
        #pragma omp parallel for num_threads(_threads)
        for(uint32_t i = 0; i < num_reads; i++)
        {
            std::set<std::string> tempset;
            #pragma omp critical
            {
                tempset.insert(umap1[pair_read[i]].begin(), umap1[pair_read[i]].end());
                tempset.insert(umap2[pair_read[i]].begin(), umap2[pair_read[i]].end());
                if (tempset.size() > 1) 
                {
                    if (seen_set.find(tempset) == seen_set.end()) 
                    {
                        seen_set.insert(tempset);
                        std::vector<std::string> tempvec{tempset.begin(), tempset.end()};
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
        FILE* ef = fopen(edgelist_file.c_str(), "w+");
        int num_sets_unitigs = vec.size();
        fprintf(stdout, "Number of unitig sets to interconnect: %zu\n", vec.size());
        for(int i = 0; i < num_sets_unitigs; i++)
        {
            std::vector<std::string> tempvec = vec[i];

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
                            if (seen_edge.find(std::make_pair(tempvec[j], tempvec[k])) == seen_edge.end())
                            {
                                #pragma omp critical
                                {
                                    fprintf(ef,"%s\t%s\n", tempvec[j].c_str(), tempvec[k].c_str());
                                    seen_edge.insert(std::make_pair(tempvec[j], tempvec[k]));
                                }
                            }
                        }
                        else
                        {
                            //check if edge already seen
                            if (seen_edge.find(std::make_pair(tempvec[k], tempvec[j])) == seen_edge.end())
                            {
                                #pragma omp critical
                                {
                                    fprintf(ef, "%s\t%s\n", tempvec[k].c_str(), tempvec[j].c_str());
                                    seen_edge.insert(std::make_pair(tempvec[k], tempvec[j]));
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

    void Kgraph::readEdgeList(const std::string& dir, const std::string& inputUnitigs)
    {
        std::string edgelist_file = dir+"/edgelist.txt";
        igraph_t graph;
        FILE* inpf = fopen(edgelist_file.c_str(), "r");
        if (inpf == nullptr) { fileNotFoundError(edgelist_file); }

        igraph_read_graph_ncol(&graph, inpf, NULL, true, IGRAPH_ADD_WEIGHTS_NO, IGRAPH_UNDIRECTED);
        fclose(inpf);
        fprintf(stdout, "GraphInfo...\n\tNumber of vertices: %d\n", (int) igraph_vcount(&graph));
        fprintf(stdout, "\tNumber of edges: %d\n", (int) igraph_ecount(&graph));

        fprintf(stdout, "Got to stage 1\n");
        
        std::map<std::string, std::string> unitigs = Kgraph::readUnitigsFile(inputUnitigs);

        fprintf(stdout, "Got to stage 2\n");

        runCore(graph, dir, unitigs);
    }

    void Kgraph::runCore(igraph_t &graph, const std::string &dir, std::map<std::string, std::string> &unitigs)
    {
        std::string kcore_file = dir+"/kcore.tsv";
        igraph_vector_int_t coreness, deg;
        igraph_vector_int_init(&coreness, 1);
        igraph_vector_int_init(&deg, 1);
        igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
        igraph_coreness(&graph, &coreness, IGRAPH_ALL);
        FILE* kcf = fopen(kcore_file.c_str(), "w+");
        int koresize = igraph_vector_int_size(&coreness);
        fprintf(kcf, "#VID\tName\tCoreness\tDegree\n");
        for(int i = 0; i < koresize; i++)
        {
            fprintf(kcf, "%d\t%s\t%d\t%d\n", i, igraph_cattribute_VAS(&graph, "name", i), (int) VECTOR(coreness)[i], (int) VECTOR(deg)[i]);
        }
        fclose(kcf);
        
        runTruss(graph, dir, coreness, koresize, unitigs); // create truss file

        igraph_vector_int_destroy(&coreness);
        igraph_vector_int_destroy(&deg);
        igraph_destroy(&graph);
    }

    void Kgraph::runTruss(igraph_t &graph, const std::string&dir, igraph_vector_int_t &coreness, int koresize, std::map<std::string, std::string> &unitigs)
    {
        igraph_vector_int_t subgraph_nodes, map, invmap;

        igraph_vector_int_init(&subgraph_nodes, 0);
        igraph_vector_int_init(&map, 0);
        igraph_vector_int_init(&invmap, 0);

        int K = (int)igraph_vector_int_max(&coreness);
        for (int i = 0; i < koresize; i++)
        {
            // get vector of nodes in maximal core
            if ((int)igraph_vector_int_get(&coreness, i) == K)
            {
                igraph_vector_int_push_back(&subgraph_nodes, i);
            }
        }
        // transform to vertex selector -- just igraph stuff...
        igraph_vs_t vids;
        igraph_vs_vector(&vids, &subgraph_nodes);

        fprintf(stdout, "BUILD K-TRUSS:\n");
        fprintf(stdout, "Selected unitigs in maximal core.\n");
        
        // get induced subgraph -- will mess up ids, so mapping required
        igraph_t subgraph;
        igraph_induced_subgraph_map(&graph, &subgraph, vids, IGRAPH_SUBGRAPH_AUTO, &map, &invmap);
        
        fprintf(stdout, "Succesfully created a K-core subgraph.\n");

        igraph_vector_int_t trussness;
        igraph_vector_int_init(&trussness, 0);
        igraph_trussness(&subgraph, &trussness);

        fprintf(stdout, "Computed trussness of edges.\n");

        int trusssize = igraph_vector_int_size(&trussness);

        std::string trussFile = dir + "/truss_unitigs.fasta";

        FILE* trussf = fopen(trussFile.c_str(), "w+");

        std::unordered_set <int> nodes;
        int threshold = K;

        while (nodes.empty())
        {
            for (int edge = 0; edge < trusssize; edge++)
            {
                if (igraph_vector_int_get(&trussness, edge) >= threshold)
                {
                    // the nodes of this edge are in the vector, add them to output

                    igraph_integer_t from, to;            
                    igraph_edge(&subgraph, edge, &from, &to);

                    nodes.insert((int)igraph_vector_int_get(&invmap, from));
                    nodes.insert((int)igraph_vector_int_get(&invmap, to));
                }
            }
            // decrement threshold, in case current truss is empty
            threshold--;
            if (threshold <= 4)
            {
                break;
            }
        }

        for (int original_node : nodes)
        {
            // print to file
            std::string unitig_header = ">Unitig_"+std::to_string(original_node);
            fprintf(trussf,"%s\n",unitig_header.c_str());
            
            auto key = unitigs.find(std::to_string(original_node));
            if (key != unitigs.end())
            {
                fprintf(trussf,"%s\n", key->second.c_str());
            }
        }
        fprintf(stdout, "%s%d%s%d%s%s\n","Found ", (int)nodes.size(), " unitigs in ", threshold+1,"-truss, saved at ", trussFile.c_str());
        
        nodes.clear();
        igraph_vector_int_destroy(&subgraph_nodes);
        igraph_vector_int_destroy(&trussness);
        igraph_vector_int_destroy(&map);
        igraph_vector_int_destroy(&invmap);
    }

    std::map<std::string, std::string> Kgraph::readUnitigsFile(const std::string& inputUnitigs)
    {
        std::map<std::string, std::string> unitigs;
        FILE* fp = fopen(inputUnitigs.c_str(),"r");
        if (fp == nullptr) { fileNotFoundError(inputUnitigs); }
        std::string unitig_num; 
        
        char* line = nullptr;
        size_t len;
        ssize_t read;
       
        while (read = getline(&line, &len, fp) != -1)
        {
            std::string uline(line); 
            if (uline.substr(0, 1) == ">")
            {
                unitig_num = uline.substr(1, uline.find(" ")-1);
                unitigs[unitig_num] = std::string("");
            } else {
                /* Need += to correctly handle multiline FASTA */
                unitigs[unitig_num] += uline.substr(0, uline.length()-1); 
            }
        }    
        return unitigs;
    }
    
    void Kgraph::combineFile(const std::string& dir, const std::string& inputUnitigs)
    {    
        std::vector<std::vector<std::string> > kcore_lines;
        std::map<std::string, std::string> unitigs;
        
        const std::string kcf = dir + "/kcore.tsv";

        FILE* f = fopen(kcf.c_str(),"r");
        if (f == nullptr) { fileNotFoundError(kcf); }
        char* line = nullptr;
        size_t len;
        ssize_t bytes_read;

        while(bytes_read = getline(&line, &len, f) != -1)
        {
            if (line[0] != '#') {
                std::vector<std::string> cur_line;
                char* token = strtok(line,"\t");
                while(token)
                {
                    std::string cur_token(token);
                    cur_line.emplace_back(cur_token);
                    token = strtok(NULL, "\t");
                }    
                kcore_lines.emplace_back(cur_line);
            }
        }
 
        unitigs = Kgraph::readUnitigsFile(inputUnitigs);
        
        const std::string cmbf = dir + "/combined.fasta";
    
        FILE* outf = fopen(cmbf.c_str(),"w+");
        for(uint64_t i = 0; i < kcore_lines.size(); i++)
        {   
            std::string header = ">Unitig_"+kcore_lines[i][1]+"|"+ kcore_lines[i][2];
            fprintf(outf, "%s\n", header.c_str());
            auto key = unitigs.find(kcore_lines[i][1]);
            if (key != unitigs.end())
            {
                fprintf(outf, "%s\n", key->second.c_str());
            }
        }
        fclose(outf);
    }
    
    void Kgraph::anomalyDetection(const std::string& dir, bool weight)
    {
        std::string edgelist_file = dir+"/edgelist.txt";                   
        igraph_t graph;                                                    
        FILE* inpf = fopen(edgelist_file.c_str(), "r");
        if (inpf == nullptr) { fileNotFoundError(edgelist_file); }                                                         
        igraph_read_graph_edgelist(&graph, inpf, 0, 0);                                                         
        igraph_lazy_adjlist_t al;                                          
        igraph_lazy_adjlist_init(&graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);                                   
        CombineCoreA::run(al, dir, weight);
    }

    double Kgraph::getMedian(std::vector<double> vec, int start, int end)
    {
       double median = 0;
       int size = end - start - 1;
       if (size % 2 == 0)
       {
           median = (vec[start + size/2 - 1] + vec[start + size/2]) / 2;
       }
       else
       {
           median = vec[start + (size - 1) / 2];
       }
       return median;

    }

    void Kgraph::splitAnomalousUnitigs(const std:: string& dir, const std:: string& inputUnitigs)
    {
        std::string anomalyFile = dir+"/top_scoring_anomalous_unitigs.txt";
        std::string backgroundFile = dir+"/low_scoring_anomalous_unitigs.txt";
        std::string coreAf =  dir+"/CoreA_anomaly.txt";
        std::map<std::string, std::string> unitigs = Kgraph::readUnitigsFile(inputUnitigs);
        std::vector<std::vector<std::string> > coreA_lines;
        FILE* inp_coreAf = fopen(coreAf.c_str(), "r");
        if (inp_coreAf == nullptr) { fileNotFoundError(coreAf); }

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
