#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <igraph.h>
#include "test.h"
#include "graph.h"
#include <chrono>
#include <ctime>
#include <sstream>
#include <tclap/CmdLine.h>
#include <utility>
#include <vector>
#include <algorithm>
#define __version__ "2.0"

int main(int argc, const char** argv)
{

    /* sync cout off */
    std::ios_base::sync_with_stdio(false);

    std::string version(__version__);
    /* check for max number of threads in the system using omp */
    int opt_threads = omp_get_max_threads();

    /* get current time info */
    time_t curr_time = time(0);
    tm* mdyt = localtime(&curr_time);
    const std::string dirname = "output_"+ std::to_string(1900 + mdyt->tm_year) +
            std::to_string(1+ mdyt->tm_mon)+ std::to_string(mdyt->tm_mday) + "_" +
            std::to_string(mdyt->tm_hour) + std::to_string(mdyt->tm_min)+
            std::to_string(mdyt->tm_sec);

    /* parse args */
    TCLAP::CmdLine cmd("KOMB: Taxonomy-oblivious characterization of metagenome "
                             "dynamics", ' ', version);
    TCLAP::ValueArg<std::string> inputArg("i", "input",
      "Input SAM file [Default: alingment1.sam]", true,
      "alignment1.sam", "string");
    TCLAP::ValueArg<std::string> input2Arg("j", "input2",
      "Second input SAM file [Default: alignment2.sam]", true,
      "alignment2.sam", "string");
   TCLAP::ValueArg<std::string> inputUnitigsArg("u", "input-unitigs",
      "FASTA file containing unitigs [Default: unitigs.fa]", true,
      "unitigs.fa", "string");
    TCLAP::ValueArg<int> readlenArg("l", "readlen",
       "Read Length (can be average) [Default: 151]", false, 151,
       "int");
    TCLAP::ValueArg<int> threadsArg("t", "threads",
       "Number of Threads [Default: Max]", false, opt_threads,
       "int");
    TCLAP::ValueArg<std::string> outputArg("o", "output",
       "Output directory [Default: output_yyyymmdd_hhmmss]", false, dirname,
       "string");
    TCLAP::SwitchArg fulgorSwitch("f", "fulgor",
       "Use Fulgor pseudoalignments instead of SAM files", cmd, false);
    cmd.add(inputArg);
    cmd.add(input2Arg);
    cmd.add(inputUnitigsArg);
    cmd.add(readlenArg);
    cmd.add(threadsArg);
    cmd.add(outputArg);

    auto begin = std::chrono::steady_clock::now(); // start timer
        
    cmd.parse(argc,argv);
    const std::string outdir = outputArg.getValue();
    const std::string input = inputArg.getValue();
    const std::string input2 = input2Arg.getValue();
    const std::string inputUnitigs = inputUnitigsArg.getValue();
    const int readlen = readlenArg.getValue();
    int threads = threadsArg.getValue();
    bool fulgor = fulgorSwitch.getValue();

    const bool isBifrost = false;
    
    /* run KOMB core */
    igraph_set_attribute_table(&igraph_cattribute_table);  // to correctly handle arbitrary unitig names
    auto kg = komb::Kgraph(threads, readlen);  // use this instantiation
    umapset umap1, umap2;
    std::unordered_map<std::string, long> unitig_id_to_vid_map;
    igraph_vector_int_t edges;
    igraph_t graph;
    long vid = 0;
    auto begin_komb = std::chrono::steady_clock::now();

//     kg.read2SAM(input, input2, umap1, unitig_id_to_vid_map, &vid, fulgor);
   //  #pragma omp parallel
   //  {
   //      #pragma omp single
   //      {
            //#pragma omp task
            kg.readSAM(input, umap1, unitig_id_to_vid_map, &vid, fulgor);
            //#pragma omp task
            kg.readSAM(input2, umap2, unitig_id_to_vid_map, &vid, fulgor);
   //      }
   //  }
    
    auto post_sam = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_sam - begin_komb).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for reading SAMs: %.3f s\n", duration);
   
    kg.getEdgeInfo(umap1, umap2);
    auto post_edgeinfo = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_edgeinfo - post_sam).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for edgeInfo: %.3f s\n", duration);
   
    unsigned vertex_count = unitig_id_to_vid_map.size();
    igraph_vector_int_init(&edges, vertex_count);  // Start at the |V| count as each vertex included has deg >= 1 
    kg.generateGraph(umap1, outdir, unitig_id_to_vid_map, edges);
    auto post_generate = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_generate - post_edgeinfo).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for generateGraph: %.3f s\n", duration);
   
    igraph_vector_int_resize_min(&edges);  // Free up any extra memory unused by edges
    kg.readEdgeList(graph, outdir, inputUnitigs, unitig_id_to_vid_map, edges);
    fprintf(stdout, "Created Kcore\n");
    auto post_core = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_core - post_generate).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for edgeInfo: %.3f s\n", duration);

//     kg.combineFile(outdir, inputUnitigs);
    auto post_combine = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_combine - post_core).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for combineFile: %.3f s\n", duration);

    kg.anomalyDetection(graph, outdir, true);
    auto post_anomaly = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            post_anomaly - post_combine).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for anomalyDetection: %.3f s\n", duration);

    fprintf(stdout,"Identified anomalous unitigs\n");
//     kg.splitAnomalousUnitigs(outdir, inputUnitigs);
    fprintf(stdout,"Created anomalouss unitigs file\n");
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - begin_komb).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for KOMB: %.3f s\n", duration);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - begin).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for analysis (sec) = %.3f \n", duration);
        
    return 0;

} // main
