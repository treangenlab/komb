//
// Created by Advait Balaji on 6/30/20.
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

typedef std::map<std::string,std::set<uint32_t> > umapset;
typedef std::vector<std::vector<uint32_t> > vvec;
typedef std::vector<std::set<uint32_t> > vecset;


namespace komb
{
    class Kgraph {
    public:
        Kgraph(int threads);
        Kgraph(int threads, int readlength);
        void readSAM(const std::string &samfile, umapset &umap);
        vvec getEdgeInfo(umapset &umap1, umapset &umap2);
        void generateGraph(std::vector<std::vector<uint32_t> > &vec,
                const std::string& dir);
        void readEdgeList(const std::string& dir);
        static void runCore(igraph_t &graph, const std::string& dir);
        void processGFA(const std::string& dir, bool weight);
        std::map<std::string, std::string> readUnitigsFile(const std::string& dir, bool isBifrost);
        void combineFile(const std::string& dir, bool isBifrost);
        void createRER(long long int& vertices, long long int& edges);
        void anomalyDetection(const std::string& dir, bool weight);
        double getMedian(std::vector<double> vec, int start, int end);
        void splitAnomalousUnitigs(const std:: string& dir, bool isBifrost);
        int _threads;
        int _readlength;
    };


    class HashSet
    {
        public:
            inline std::size_t operator()(const std::set<uint32_t> &s) const
            {
                std::vector<int> sums(s.size(), 0);
                std::vector<int> products(s.size());
                std::partial_sum(s.begin(), s.end(), sums.begin());
                std::partial_sum(s.begin(), s.end(), products.begin(),
                             std::multiplies<int>());
                return sums[s.size() - 1] + products[s.size() - 1];
            }
    };

    class HashPairSet
    {
        public:
            inline std::size_t  operator()(const std::pair<uint32_t, uint32_t> &p ) const
            {
                return p.first + p.second + (p.first * p.second);
            }
    };

}

typedef std::set<std::set<uint32_t> > usset;
typedef std::set<std::pair<uint32_t, uint32_t> > uspair;
#endif //KOMB_GRAPH_H
