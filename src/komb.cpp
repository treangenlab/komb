#include <cstdlib>
#include <iostream>
#include <CompactedDBG.hpp>
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
#define __version__ "1.0"


CDBG_Build_opt setOpts(int numthreads, bool deleteisolated, bool outputgfa, bool cliptips,
        bool verbose)
{
    CDBG_Build_opt opt;
    opt.nb_threads = numthreads;
    opt.deleteIsolated = deleteisolated;
    opt.outputGFA = outputgfa;
    opt.clipTips = cliptips;
    opt.verbose = verbose;

    return opt;
}

void runSpades(std::vector<std::string>& reads, const std::string& outdir, int& threads)
{
    int return_val;
    std::string readfile1 = "\t \""+reads[0]+"\",\n";
    std::string readfile2 = "\t \""+reads[1]+"\",\n";
    std::string yml_file = outdir + "/input.yaml";
    fprintf(stdout, "Writing the Spades yaml file\n");
    FILE* f = fopen(yml_file.c_str(),"w+");
    fprintf(f,"[\n");
    fprintf(f,"\t{\n");
    fprintf(f,"\torientation: \"fr\",\n");
    fprintf(f,"\ttype: \"paired-end\",\n");
    fprintf(f,"\tright reads: [\n");
    fprintf(f,"%s",readfile1.c_str());
    fprintf(f,"\t],\n");
    fprintf(f,"\tleft reads: [\n");
    fprintf(f,"%s",readfile2.c_str());
    fprintf(f,"\t],\n");
    fprintf(f,"\t}\n");
    fprintf(f,"]");

    fprintf(stdout, "Running spades\n");
    std::string spades_run = "spades-gbuilder "+ outdir+"/input.yaml " +
            outdir+"/output.gfa -k 31 -t "+ std::to_string(threads) + " --gfa";
    return_val = std::system(spades_run.c_str());
}


void runBowtie2(const int& numhits, const int& readlen, std::vector<std::string>& reads,
        const std::string& outdir, const std::string& fast_x, const int& kmersize, int& threads)
{
    int return_val;
    return_val = std::system("bowtie2-build unitigs.fasta idx");
    fprintf(stdout,"Built Bowtie2 index for unitigs\n");
    std::string align1 = "bowtie2 " + fast_x + " -x idx --sensitive -k " + std::to_string(numhits) +
            " -3 " + std::to_string(readlen-kmersize+1) + " -U " + reads[0] +
            " -p " + std::to_string(threads) + " --no-unal > alignment1.sam";
    std::string align2 = "bowtie2 " + fast_x + " -x idx --sensitive -k " + std::to_string(numhits) +
                         " -3 " + std::to_string(readlen-kmersize+1) + " -U " + reads[1] +
                         " -p " + std::to_string(threads) + " --no-unal > alignment2.sam";
    return_val = std::system(align1.c_str());
    fprintf(stdout,"Aligned first set of reads\n");
    return_val = std::system(align2.c_str());
    fprintf(stdout,"Aligned second set of reads\n");
    return_val = std::system("rm -fr idx*");
    fprintf(stdout,"Cleaned index files and completed Bowtie2\n");
    std::string move_unitigs = "mv unitigs.fasta "+outdir;
    std::string move_alignments = "mv *.sam "+outdir;
    return_val = std::system(move_unitigs.c_str());
    return_val = std::system(move_alignments.c_str());
}

inline bool checkGC(std::string& unitig)
{
        double gc_content = (std::count(unitig.begin(), unitig.end(), 'G')+std::count(unitig.begin(), unitig.end(), 'C'))/unitig.length();
        if (0.1 < gc_content < 0.9)
        {
           return true;
        }
        return false;
}

void fixAbyssUnitigs()
{
     std::vector<std::pair<std::string, std::string> > unitigs;         
     const std::string utgf = "temp-unitigs.fa";                          
                                                                      
     FILE* fp = fopen(utgf.c_str(),"r");                                
     if (fp == nullptr) { exit(EXIT_FAILURE);}                          
     std::string unitig_num;                                            
     char* line = nullptr;                                              
     size_t len;                                                        
     ssize_t read;                                                      
     while (read = getline(&line,&len,fp) != -1)                        
     {                                                                  
       std::string uline(line);                                         
       if (uline.substr(0,1) == ">")                                    
       {                                                                
          unitig_num = uline.substr(uline.find(" ")+1);//,uline.find("\n")-5); // 5 just seems to work here!
          unitig_num.erase(std::remove(unitig_num.begin(), unitig_num.end(), '\n'), unitig_num.end());
       }                                                                
       else                                                             
       {  
          std::string useq = uline.substr(0,uline.length()-1);
          if (checkGC(useq))
          {
            unitigs.emplace_back(std::make_pair(unitig_num, useq));      
          }
       }                                                                
     }                                                                  
                                                                         
    FILE* wp = fopen("unitigs.fasta","w+");                        
    for (int i = 0; i < unitigs.size(); i++)                           
    {                                                                  
      fprintf(wp,">%d %s\n",i,unitigs[i].first.c_str());             
      fprintf(wp,"%s\n",unitigs[i].second.c_str());                  
    }                                                                  
    fclose(wp);          
}

void runAbyss(std::vector<std::string>& reads, const int& kmersize, int& processes)
{
    int return_val;
    std::string runabyss = "abyss-pe np=" + std::to_string(processes) + " name=temp k=" + std::to_string(kmersize) +
            " in='" + reads[0] + " " + reads[1] + "' unitigs";
    return_val = std::system(runabyss.c_str());
    fixAbyssUnitigs();
    //std::string rename_unitigs = "cp temp-unitigs.fa unitigs.fasta";
    //return_val = std::system(rename_unitigs.c_str());
    return_val = std::system("rm -rf temp-* coverage.hist");      
    fprintf(stdout, "Cleaned Abyss temp, addditional and hist files\n");      
}

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
    TCLAP::ValueArg<std::string> readArg("r","reads",
      "Paired-read file separated by ',' [Default: read1.fq,read2.fq]",true,
      "read1.fq,read2.fq","string");
    TCLAP::ValueArg<int> readlenArg("l","readlen",
       "Read Length (can be average) [Default: 100]",false,100,
       "int");
    TCLAP::ValueArg<int> numhitsArg("n","numhits",
        "Bowtie2 maximum number of hits [Default: 1000]",false,1000,
        "int");
    TCLAP::ValueArg<int> kmerArg("k","kmer",
        "Kmer size for Abyss, Bifrost uses 31 [Default: 31]",false,31,
        "int");
    TCLAP::ValueArg<int> threadsArg("t","threads",
       "Number of Threads [Default: Max]",false,opt_threads,
       "int");
    TCLAP::SwitchArg spadesSwitch("s","spades",
       "Runs spades and uses GFA graph instead of bifrost + bowtie2 [Default: run abyss]",cmd, false);
    //switch implicitly adds cmd
    TCLAP::SwitchArg bifrostSwitch("b","bifrost",
       "Run bifrost instead of abyss [Default: run abyss]",cmd, false);
    //switch implicitly adds cmd
    TCLAP::SwitchArg fastaSwitch("f","fasta",
       "Reads provided are fasta files [Default: fastq]",cmd,false);
    //switch implicitly adds cmd
    TCLAP::SwitchArg alignfileSwitch("a","alignment",
       "Keep alignment files [Default: delete alignment]",cmd,false);
    //switch implicitly adds cmd
    TCLAP::ValueArg<std::string> outputArg("o","output",
       "Output directory [Default: output_yyyymmdd_hhmmss]", false,dirname,
       "string");
    cmd.add(readArg);
    cmd.add(readlenArg);
    cmd.add(numhitsArg);
    cmd.add(kmerArg);
    cmd.add(threadsArg);
    cmd.add(outputArg);

    auto begin = std::chrono::steady_clock::now(); // start timer
        
    cmd.parse(argc,argv);
    std::string readargs = readArg.getValue();
    const int readlen = readlenArg.getValue();
    const int kmer = kmerArg.getValue();
    int threads = threadsArg.getValue();
    const int numhits = numhitsArg.getValue();
    const std::string outdir = outputArg.getValue();
    bool spades = spadesSwitch.getValue();
    bool isBifrost = bifrostSwitch.getValue();
    bool isFasta = fastaSwitch.getValue();
    bool keepAlign = alignfileSwitch.getValue();

    /* get the readfiles */
    std::vector<std::string> reads;
    std::stringstream ss(readargs);
    while( ss.good() )
    {
        string substr;
        getline( ss, substr, ',' );
        reads.emplace_back( substr );
    }

    /* sanity checks for arguments readfile*/
    if (reads.size() != 2)
    {
        fprintf(stderr, "Reads not provided/recognized or greater than 2 files\n");
        exit(1);
    }
        
    /* sanity check if kmer size is greater thatn the read length */
    if (kmer > readlen)
    {
        fprintf(stderr, "Kmer size is greater than read length\n");
        exit(1);
    }
        
    /* raise warning if kmersize is greater than 75% of read length */
    if (kmer >= 0.75*readlen)
    {
        fprintf(stderr, "[Warning] Kmer size is greater than 75%s of the read length\n","%");
    }
        
    /* create output dir */
    const std::string create_dir = "mkdir " + outdir;
    int return_val = std::system(create_dir.c_str());

    if (spades)
    {
        runSpades(reads, outdir, threads);
        auto kg = komb::Kgraph(threads, readlen);
        //int weight = 10; /* fix weight as 10 for now */
        kg.processGFA(outdir, true);
        exit(0);
    }
        
    /* set kmer size */
    const int kmersize = isBifrost ? 31 : kmer;

    if (!isBifrost)
    {
        runAbyss(reads, kmersize, threads);
        fprintf(stdout, "Created unitigs from Abyss\n");
    }
    else
    {
        /* instantiate and run bifrost */
        
        CompactedDBG<> cdbg; // cdbg init
        assert(cdbg.getK() == 31); // assert that k-mer size is default for now
        CDBG_Build_opt opt = setOpts(threads, true, false, true,
                                 true);
        for (auto &read : reads)
        {
                opt.filename_seq_in.emplace_back(read);
        }
        opt.prefixFilenameOut = "unitigs";

        cdbg.build(opt); //build the DBG
        cdbg.simplify(opt.deleteIsolated, opt.clipTips,
                  opt.verbose); //delete false unitigs
        cdbg.write(opt.prefixFilenameOut, opt.nb_threads, 
                   opt.outputGFA, opt.verbose);
        fprintf(stdout, "Created unitigs from bifrost\n");
    }
     
    /* set fasta or fastq for bowtie */
    const std::string fast_x = isFasta ? "-f" : "-q";

    /* run bowtie2 */
    runBowtie2(numhits, readlen, reads, outdir, fast_x, kmersize, threads);

    /* run KOMB core */

    auto kg = komb::Kgraph(threads, readlen); //use this instantiation

    const std::string alignment1 = outdir + "/alignment1.sam";
    const std::string alignment2 = outdir + "/alignment2.sam";
    umapset umap1, umap2;

    #pragma omp parallel
    {
        #pragma omp single
        {
            //#pragma omp task
            kg.readSAM(alignment1, umap1);
            //#pragma omp task
            kg.readSAM(alignment2, umap2);
        }

    }

    fprintf(stdout, "Read SAM files\n");
    if (!keepAlign) //remove alignment files as they take up too much space, unless the user wants it!
    {
        std::string remove_sam = "rm -rf "+outdir+"/*.sam";
        int return_val = std::system(remove_sam.c_str());
    }
    vvec edgeinfo = kg.getEdgeInfo(umap1, umap2);
    fprintf(stdout, "Processed edgeinfo\n");
    kg.generateGraph(edgeinfo, outdir);
    fprintf(stdout, "Created edgelist\n");
    kg.readEdgeList(outdir);
    fprintf(stdout, "Created Kcore\n");
    kg.combineFile(outdir, isBifrost);
    //int weight = 10; /* fix weight as 10 for now */
    kg.anomalyDetection(outdir, true);
    fprintf(stdout,"Identified anomalous unitigs\n");
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - begin).count() / 1000000.0;
    fprintf(stdout, "\nTime elapsed for analysis (sec) = %.2f \n", duration);
        
    return 0;

} // main
