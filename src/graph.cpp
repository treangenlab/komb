//
// Created by Advait Balaji on 06/30/2020.
// Modified by Nicolae Sapoval on 06/25/2023.
// Modified by Marko Tanevski on 07/07/2023
//

#include <cstring>
#include "graph.h"
#include <omp.h>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <functional>
#include <utility>

using namespace gfa;

namespace komb
{
    struct hash_pair final {
        template<class TFirst, class TSecond>
        size_t operator()(const std::pair<TFirst, TSecond>& p) const noexcept {
            uintmax_t hash = std::hash<TFirst>{}(p.first);
            hash <<= sizeof(uintmax_t) * 4;
            hash ^= std::hash<TSecond>{}(p.second);
            return std::hash<uintmax_t>{}(hash);
        }
    };

    template <typename T>
    bool is_subset_of(const std::unordered_set<T>& a, const std::unordered_set<T>& b)
    {
        // return true if all members of a are also in b
        if (a.size() > b.size())
            return false;

        auto const not_found = b.end();
        for (auto const& element: a)
            if (b.find(element) == not_found)
                return false;

        return true;
    }

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

    
    void Kgraph::readSAM(const std::string &samfile, umapset &umap, std::unordered_map<std::string, long> &unitig_id_to_vid_map, long *vid, bool fulgor)
    {
        FILE* f = fopen(samfile.c_str(), "r");
        if (f == nullptr) { fileNotFoundError(samfile); }
        char* file_buffer;
        size_t len;
        uint64_t file_size;
        uint64_t bytes_read;

        // Get file size
        fseek(f, 0, SEEK_END);
        file_size = ftell(f);
        rewind(f);

        // Allocate memory and slurp the whole file in
        file_buffer = (char *) malloc(sizeof(char)*file_size+1);
        if (file_buffer == NULL) 
        { 
            std::cerr << "Could not allocate memory for reading " << samfile << std::endl;
            exit(EXIT_FAILURE);
        }
        bytes_read = fread(file_buffer, 1, file_size, f);
        if (bytes_read != file_size) 
        {
            std::cerr << "Encountered error while reading " << samfile << std::endl;
            exit(EXIT_FAILURE);
        }

        // Process file in parallel (based on: https://stackoverflow.com/questions/16812302/openmp-while-loop-for-text-file-reading-and-using-a-pipeline)
        std::vector<uint64_t> *position;
        umapset local_umaps[_threads];
        #pragma omp parallel num_threads(_threads)
        {
            const int ithread = omp_get_thread_num();  // Get thread number to index into local umaps
            #pragma omp single 
            {
                position = new std::vector<uint64_t>[_threads];
                position[0].push_back(0);
            }

            #pragma omp for 
            for(uint64_t i=0; i<bytes_read; ++i) 
            {
                if(file_buffer[i] == '\n' || file_buffer[i] == '\0') 
                {
                    position[ithread].push_back(i);
                }
            }

            for (uint64_t i=1; i<position[ithread].size(); ++i)
            {
                // position points to location of '\n' so we need to offset by 1, unless it is the start aka 0
                uint64_t start_pos = position[ithread][i-1]+1;
                if(start_pos == 1) start_pos = 0;
                char *line = &file_buffer[start_pos];
                if(line[0] != '@')
                {
                    char* saveptr;
                    char* token = strtok_r(line, "\t", &saveptr);
                    std::string read(token);
                    uint8_t token_count = 0;
                    while(token && (token_count < 2))
                    {
                        token = strtok_r(NULL, "\t", &saveptr);
                        token_count++;
                    }
                    std::string unitig(token);
                    if (unitig.compare("*") != 0) 
                    {  // '*' indicates an unmapped read in SAM
                        local_umaps[ithread][read.substr(1, read.find('/'))].insert(unitig);
                    }
                }
            }
        }
        
        // Now reduce per thread umaps into a single one
        for (int i=0; i<_threads; ++i)
        {
            for (auto &it : local_umaps[i]) 
            {
                for (auto unitig_it = it.second.begin(); unitig_it!=it.second.end(); ++unitig_it)
                {
                    umap[it.first].insert((*unitig_it));
                    if (auto uid2vid_it = unitig_id_to_vid_map.find((*unitig_it)); 
                        uid2vid_it == unitig_id_to_vid_map.end()) 
                    {
                        unitig_id_to_vid_map[(*unitig_it)] = (*vid)++;  
                    }
                }
            }
        }
    }

    void Kgraph::getEdgeInfo(umapset &umap1, umapset &umap2)
    {
        /* Merge paired-end read information into umap1:
            each key (read name) in umap1 will contain all unitig_IDs
            that are connected by that read or via paired-end information.
            This is OK since in a HUG that would always form a clique.
         */
        unsigned n = umap1.bucket_count();  
        for(unsigned i=0; i<n; ++i)
        {
            for(auto it1=umap1.begin(i); it1!=umap1.end(i); ++it1)
            {
                std::string read = (*it1).first;
                if (auto it2 = umap2.find(read); it2 != umap2.end())
                {
                    umap1[read].merge((*it2).second);  // Merge will skip unitig names that are already in umap1[read]
                    umap2.erase(it2);  // All information now extracted, delete the element (reduces size of umap2 for the next step)
                }
            }
        }
        
        // Now add any cliques that arise from singleton mappings in umap2
        umap1.insert(umap2.begin(), umap2.end());
        
        // Clean up after ourselves; .swap() is currently preferred method
        umapset().swap(umap2);
    }

    void Kgraph::generateGraph(umapset& umap, const std::string& dir, 
                               std::unordered_map<std::string, long> &unitig_id_to_vid_map,
                               igraph_vector_int_t &edges)
    {
        // std::string edgelist_file = dir+"/edgelist.txt";
        // FILE* ef = fopen(edgelist_file.c_str(), "w+");
        uspair seen_edges;
        auto begin_vec = std::chrono::steady_clock::now();

        std::vector<std::vector<std::string>> cliques;
        cliques.reserve(umap.size());
        for(auto &it : umap)
        {
            std::vector<std::string> clique_nodes(it.second.begin(), it.second.end());
            cliques.push_back(clique_nodes);
        }
        umapset().swap(umap);

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - begin_vec).count() / 1000000.0;
        fprintf(stdout, "\nTime elapsed for converting umapset to vec<vec>: %.3f s\n", duration);

        auto begin_insert = std::chrono::steady_clock::now();

        std::vector<long> local_edges[_threads];
        std::unordered_set<std::pair<long, long>, hash_pair> local_seen_edges[_threads];

        #pragma omp parallel num_threads(_threads) shared(unitig_id_to_vid_map)
        {
            const int ithread = omp_get_thread_num();  // Get thread number to index into local edge vectors
            #pragma omp for 
            for(unsigned clique_id=0; clique_id<cliques.size(); ++clique_id)  
            {
                for(unsigned i=0; i<cliques[clique_id].size(); ++i) 
                {
                    for(unsigned j=i+1; j<cliques[clique_id].size(); ++j)
                    {
                        std::string uid1 = cliques[clique_id][i];
                        std::string uid2 = cliques[clique_id][j];
                        std::pair<long, long> edge = std::make_pair(unitig_id_to_vid_map[uid1], unitig_id_to_vid_map[uid2]);

                        if (local_seen_edges[ithread].find(edge) == local_seen_edges[ithread].end()) 
                        {
                            local_edges[ithread].push_back(unitig_id_to_vid_map[uid1]);
                            local_edges[ithread].push_back(unitig_id_to_vid_map[uid2]);
                            local_seen_edges[ithread].insert(edge);
                        }
                    }
                }
            } 
        }

        #pragma omp parallel for num_threads(_threads)
        for(int i=0; i<_threads; ++i) 
        {
            std::unordered_set<std::pair<long, long>, hash_pair>().swap(local_seen_edges[i]);
        }

        duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - begin_insert).count() / 1000000.0;
        fprintf(stdout, "\nTime elapsed for constructing local edges: %.3f s\n", duration);

        uint64_t total_edge_countx2 = 0;
        uint64_t edge_offsets[_threads];
        for(int i=0; i<_threads; ++i) 
        {
            edge_offsets[i] = total_edge_countx2;
            total_edge_countx2 += local_edges[i].size();
        } 
        igraph_vector_int_resize(&edges, total_edge_countx2); 
        
        #pragma omp parallel for shared(edges)
        for(int i=0; i<_threads; ++i) 
        {
            for (int j=0; j<local_edges[i].size(); j+=2)
            {
                igraph_vector_int_set(&edges, edge_offsets[i]+j, local_edges[i][j]);
                igraph_vector_int_set(&edges, edge_offsets[i]+j+1, local_edges[i][j+1]);
            }
        }

        // fclose(ef);
        uspair().swap(seen_edges);
    }

    void Kgraph::readEdgeList(igraph_t &graph, const std::string& dir, const std::string& inputUnitigs, 
                              std::unordered_map<std::string, long> &unitig_id_to_vid_map,
                              igraph_vector_int_t &edges)
    {
        // std::string edgelist_file = dir+"/edgelist.txt";
        // igraph_t graph;
        // igraph_strvector_t vid_to_uid;
        // FILE* inpf = fopen(edgelist_file.c_str(), "w");
        // if (inpf == nullptr) { fileNotFoundError(edgelist_file); }

        std::vector<std::string> vid_to_uid (unitig_id_to_vid_map.size());
        for (auto & uid2vid : unitig_id_to_vid_map) 
        {
            vid_to_uid[uid2vid.second] = uid2vid.first;
        }

        auto begin_graph = std::chrono::steady_clock::now();
        // #pragma omp parallel num_threads(2)
        // {
        //     #pragma omp sections
        //     {
        //         #pragma omp section
        //         {
        //             long num_vertices = unitig_id_to_vid_map.size();
        //             igraph_strvector_init(&vid_to_uid, num_vertices);
        //             for (auto & uid2vid : unitig_id_to_vid_map) 
        //             {
        //                 igraph_strvector_set_len(&vid_to_uid, uid2vid.second, uid2vid.first.c_str(), uid2vid.first.length());
        //             }
                    igraph_create(&graph, &edges, unitig_id_to_vid_map.size(), 0 /* undirected */);
                //     SETVASV(&graph, "name", &vid_to_uid);
                // }
        //         #pragma omp section
        //         {
        //             for(int i=0; i<igraph_vector_int_size(&edges); i+=2) 
        //             {
        //                 fprintf(inpf, "%ld\t%ld\n", VECTOR(edges)[i], VECTOR(edges)[i+1]);
        //             }
        //         }
        //     }
        // }
        
        // igraph_read_graph_ncol(&graph, inpf, NULL, true, IGRAPH_ADD_WEIGHTS_NO, IGRAPH_UNDIRECTED);
        // fclose(inpf);
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - begin_graph).count() / 1000000.0;
        fprintf(stdout, "\nTime elapsed for initializing igraph graph: %.3f s\n", duration);

        auto begin_simplify = std::chrono::steady_clock::now();
        igraph_simplify(&graph, /*multiple=*/ true, /*loops=*/ true, /*edge_comb=*/ NULL); // discard multiple edges and self loops
        duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - begin_simplify).count() / 1000000.0;
        fprintf(stdout, "\nTime elapsed for simplifying graph: %.3f s\n", duration);
        
        fprintf(stdout, "GraphInfo...\n\tNumber of vertices: %d\n", (int) igraph_vcount(&graph));
        fprintf(stdout, "\tNumber of edges: %d\n", (int) igraph_ecount(&graph));
        
        std::unordered_map<std::string, std::string> unitigs = Kgraph::readUnitigsFile(inputUnitigs);

        auto begin_kcore = std::chrono::steady_clock::now();
        runCore(graph, dir, unitigs, vid_to_uid);
        duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - begin_kcore).count() / 1000000.0;
        fprintf(stdout, "\nTime elapsed doing K-core decomposition: %.3f s\n", duration);
    }

    void Kgraph::runCore(igraph_t &graph, const std::string &dir, std::unordered_map<std::string, std::string> &unitigs,
                         std::vector<std::string> vid_to_uid_map)
    {
        std::string kcore_file = dir+"/kcore.tsv";
        igraph_vector_int_t coreness, deg, subgraph_nodes;
        igraph_vector_int_init(&subgraph_nodes, 0);
        igraph_vector_int_init(&coreness, 1);
        igraph_vector_int_init(&deg, 1);
        igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
        igraph_coreness(&graph, &coreness, IGRAPH_ALL);
        FILE* kcf = fopen(kcore_file.c_str(), "w+");
        int koresize = igraph_vector_int_size(&coreness);
        const int max_coreness = (int)igraph_vector_int_max(&coreness);
        fprintf(kcf, "#VID\tName\tCoreness\tDegree\n");
        for(int i = 0; i < koresize; i++)
        {
            if ((int) VECTOR(coreness)[i] == max_coreness) 
            {
                igraph_vector_int_push_back(&subgraph_nodes, i);
            }
            fprintf(kcf, "%d\t%s\t%d\t%d\n", i, vid_to_uid_map[i].c_str(), (int) VECTOR(coreness)[i], (int) VECTOR(deg)[i]);
        }
        fclose(kcf);
        
        // runTruss(graph, dir, subgraph_nodes, max_coreness, unitigs); // create truss file

        igraph_vector_int_destroy(&subgraph_nodes);
        igraph_vector_int_destroy(&coreness);
        igraph_vector_int_destroy(&deg);
        // igraph_destroy(&graph);
    }

    void Kgraph::runTruss(igraph_t &graph, const std::string&dir, igraph_vector_int_t &subgraph_nodes, const int K, std::unordered_map<std::string, std::string> &unitigs)
    {
        igraph_vector_int_t map, invmap;

        igraph_vector_int_init(&map, 0);
        igraph_vector_int_init(&invmap, 0);

        // transform to vertex selector -- just igraph stuff...
        igraph_vs_t vids;
        igraph_vs_vector(&vids, &subgraph_nodes);

        fprintf(stdout, "BUILDING K-TRUSS:\n");
        fprintf(stdout, "Selected unitigs in maximal core.\n");
        
        // get induced subgraph -- will mess up ids, so mapping required
        igraph_t subgraph;
        igraph_induced_subgraph_map(&graph, &subgraph, vids, IGRAPH_SUBGRAPH_AUTO, &map, &invmap);
        
        fprintf(stdout, "%s%d%s%d%s\n", "Succesfully created a ", K, "-core subgraph, with ", (int) igraph_ecount(&subgraph), " edges.");

        igraph_vector_int_t trussness;
        igraph_vector_int_init(&trussness, 0);
        igraph_trussness(&subgraph, &trussness);

        fprintf(stdout, "Computed trussness of edges.\n");

        int trusssize = igraph_vector_int_size(&trussness);

        std::string trussFile = dir + "/truss_unitigs.fasta";

        FILE* trussf = fopen(trussFile.c_str(), "w+");

        std::unordered_set <int> nodes;
        int threshold = (int)igraph_vector_int_max(&trussness);

        // while (nodes.empty())
        // {
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
        //     // decrement threshold, in case current truss is empty
        //     threshold--;
        //     if (threshold <= 4)
        //     {
        //         break;
        //     }
        // }

        // print to file
        for (int original_node : nodes)
        {
            // Looks up the original unitig name based on the vertex ID
            std::string unitig_name = igraph_cattribute_VAS(&graph, "name", original_node); 
            std::string unitig_header = ">Unitig_"+unitig_name;
            fprintf(trussf,"%s\n",unitig_header.c_str());
            
            auto key = unitigs.find(unitig_name);
            if (key != unitigs.end())
            {
                fprintf(trussf,"%s\n", key->second.c_str());
            }
        }
        fprintf(stdout, "%s%d%s%d%s%s\n","Found ", (int)nodes.size(), " unitigs in ", threshold+1,"-truss, saved at ", trussFile.c_str());
        
        nodes.clear();
        igraph_vector_int_destroy(&trussness);
        igraph_vector_int_destroy(&map);
        igraph_vector_int_destroy(&invmap);
    }

    std::unordered_map<std::string, std::string> Kgraph::readUnitigsFile(const std::string& inputUnitigs)
    {
        std::unordered_map<std::string, std::string> unitigs;
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
        std::unordered_map<std::string, std::string> unitigs;
        
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
    
    void Kgraph::anomalyDetection(igraph_t &graph, const std::string& dir, bool weight)
    {
        // std::string edgelist_file = dir+"/edgelist.txt";                   
        // igraph_t graph;                                                    
        // FILE* inpf = fopen(edgelist_file.c_str(), "r");
        // if (inpf == nullptr) { fileNotFoundError(edgelist_file); }                                                         
        // igraph_read_graph_edgelist(&graph, inpf, 0, 0);                                                         
        igraph_lazy_adjlist_t al;                                          
        igraph_lazy_adjlist_init(&graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);                                   
        CombineCoreA::run(al, dir, weight);
        igraph_destroy(&graph);
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
        std::unordered_map<std::string, std::string> unitigs = Kgraph::readUnitigsFile(inputUnitigs);
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
