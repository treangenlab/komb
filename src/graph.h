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

typedef std::map<std::string, std::set<std::string>> umapset;
typedef std::vector<std::vector<std::string> > vvec;
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
        void readSAM(const std::string &samfile, umapset &umap, bool fulgor);
        vvec getEdgeInfo(umapset &umap1, umapset &umap2);
        void generateGraph(vvec &vec, const std::string& dir);
        void readEdgeList(const std::string &dir, const std::string &inputUnitigs);
        static void runCore(igraph_t &graph, const std::string& dir, std::map<std::string, std::string> &unitigs);
        static void runTruss(igraph_t &graph, const std::string &dir, igraph_vector_int_t &subgraph_nodes, const int K, std::map<std::string, std::string> &unitigs);
        std::map<std::string, std::string> readUnitigsFile(const std::string& inputUnitigs);
        void combineFile(const std::string& dir, const std::string& inputUnitigs);
        void anomalyDetection(const std::string& dir, bool weight);
        double getMedian(std::vector<double> vec, int start, int end);
        void splitAnomalousUnitigs(const std:: string& dir, const std:: string& inputUnitigs);
    };
}

typedef std::set<std::set<std::string>> usset;
typedef std::set<std::pair<std::string, std::string>> uspair;
#endif //KOMB_GRAPH_H
