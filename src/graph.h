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
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <igraph.h>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "gfa.h"
#include "CombineCoreA.h"

typedef std::unordered_map<std::string, std::unordered_set<std::string>> umapset;
typedef std::vector<std::vector<std::string>> vvec;
typedef std::vector<std::set<uint32_t>> vecset;
typedef std::set<std::unordered_set<std::string>> usset;  // Have to use set, since there is no easy default HASH(unordered_set<string>);

struct StringPairHash
{
    auto operator()(const std::pair<std::string, std::string> &p) const -> size_t 
    {
        return std::hash<std::string>{}(p.first + p.second);
    }
};
typedef std::unordered_set<std::pair<std::string, std::string>, StringPairHash> uspair;

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
        void readSAM(const std::string &samfile, umapset &umap, std::unordered_map<std::string, long> &unitig_id_to_vid_map, long *vid, bool fulgor);
        void getEdgeInfo(umapset &umap1, umapset &umap2);
        void generateGraph(umapset &umap, const std::string& dir, std::unordered_map<std::string, long> &unitig_id_to_vid_map, igraph_vector_int_t &edges);
        void readEdgeList(igraph_t &graph, const std::string &dir, const std::string &inputUnitigs, std::unordered_map<std::string, long> &unitig_id_to_vid_map, igraph_vector_int_t &edges);
        static void runCore(igraph_t &graph, const std::string& dir, std::unordered_map<std::string, std::string> &unitigs, std::vector<std::string> vid_to_uid);
        static void runTruss(igraph_t &graph, const std::string &dir, igraph_vector_int_t &subgraph_nodes, const int K, std::unordered_map<std::string, std::string> &unitigs);
        std::unordered_map<std::string, std::string> readUnitigsFile(const std::string& inputUnitigs);
        void combineFile(const std::string& dir, const std::string& inputUnitigs);
        void anomalyDetection(igraph_t &graph, const std::string& dir, bool weight);
        double getMedian(std::vector<double> vec, int start, int end);
        void splitAnomalousUnitigs(const std:: string& dir, const std:: string& inputUnitigs);
    };
}
#endif //KOMB_GRAPH_H
