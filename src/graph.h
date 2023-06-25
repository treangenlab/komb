//
// Created by Advait Balaji on 06/30/2020.
// Modified by Nicolae Sapoval on 06/25/2023.
//

#ifndef KOMB_GRAPH_H
#define KOMB_GRAPH_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdio>
#include <set>
#include <map>
#include <vector>
#include <igraph.h>
#include <algorithm>
#include <numeric>
#include "gfa.h"
#include "CombineCoreA.h"

typedef std::map<std::string,std::set<uint64_t> > umapset;
typedef std::vector<std::vector<uint32_t> > vvec;
typedef std::vector<std::set<uint32_t> > vecset;

namespace komb
{
    class Kgraph {
        uint32_t _threads;
        uint64_t _readlength;

    private:
        void fileNotFoundError(const std::string& path);

    public:
        Kgraph(uint32_t threads);
        Kgraph(uint32_t threads, uint64_t readlength);
        void readSAM(const std::string &samfile, umapset &umap);
        vvec getEdgeInfo(umapset &umap1, umapset &umap2);
        void generateGraph(std::vector<std::vector<uint32_t> > &vec, const std::string& dir);
        void readEdgeList(const std::string& dir);
        static void runCore(igraph_t &graph, const std::string& dir);
        void processGFA(const std::string& dir, bool weight);
        std::map<std::string, std::string> readUnitigsFile(const std::string& inputUnitigs);
        void combineFile(const std::string& dir, const std::string& inputUnitigs);
        void createRER(long long int& vertices, long long int& edges);
        void anomalyDetection(const std::string& dir, bool weight);
        double getMedian(std::vector<double> vec, int start, int end);
        void splitAnomalousUnitigs(const std:: string& dir, const std:: string& inputUnitigs);
    };
}

typedef std::set<std::set<uint32_t> > usset;
typedef std::set<std::pair<uint32_t, uint32_t> > uspair;
#endif //KOMB_GRAPH_H
