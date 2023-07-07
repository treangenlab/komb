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

   /* create output dir */
   //  const std::string create_dir = "mkdir " + outdir;
   //  int return_val = std::system(create_dir.c_str());
   //  if (return_val != 0) { 
   //    std::cerr << "Could not create output directory: " << create_dir << ". Exiting..." << std::endl;
   //    exit(EXIT_FAILURE);
   //  }

    /* run KOMB core */
    igraph_set_attribute_table(&igraph_cattribute_table);  // to correctly handle arbitrary unitig names
    auto kg = komb::Kgraph(threads, readlen);  // use this instantiation
    umapset umap1, umap2;

    #pragma omp parallel
    {
        #pragma omp single
        {
            //#pragma omp task
            kg.readSAM(input, umap1, fulgor);
            //#pragma omp task
            kg.readSAM(input2, umap2, fulgor);
        }
    }

   //  fprintf(stdout, "Read SAM files\n");
    
    auto begin_komb = std::chrono::steady_clock::now();
    vvec edgeinfo = kg.getEdgeInfo(umap1, umap2);
   //  fprintf(stdout, "Processed edgeinfo\n");
    kg.generateGraph(edgeinfo, outdir);
   //  fprintf(stdout, "Created edgelist\n");
    kg.readEdgeList(outdir);
    fprintf(stdout, "Created Kcore\n");
    kg.combineFile(outdir, inputUnitigs);
    //int weight = 10; /* fix weight as 10 for now */
    kg.anomalyDetection(outdir, true);
    fprintf(stdout,"Identified anomalous unitigs\n");
    kg.splitAnomalousUnitigs(outdir, inputUnitigs);
    fprintf(stdout,"Created anomalouss unitigs file\n");
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - begin_komb).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for KOMB: %.3f ms\n", duration);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - begin).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for analysis (msec) = %.3f \n", duration);
        
    return 0;

} // main
