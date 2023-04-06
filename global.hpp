#ifndef _GLOBAL_HPP_
#define _GLOBAL_HPP_
#include "SSA/hypergraph.hpp"
#include "GA_src/rng.h"

Graph g_g;
HyperGraph g_hg;

// fitness calculation
std::vector<int> g_nodeDegree, // record degree of node
    g_nodeDegree_cpy;          // for calculating the fitness
std::vector<bool> g_edgeMark;  //

size_t num_refresh = 10000;
// repair function
std::vector<std::vector<int>> g_edgeList3;

inline void refresh(bool lt)
{
    g_hg.clearEdges();
    addHyperedge(g_g, g_hg, 1, num_refresh, lt);

    for (vector<int>::const_iterator it = g_hg.getNodeMapVec().cbegin(); it != g_hg.getNodeMapVec().cend(); it++)
        g_nodeDegree[*it - 1] = g_hg.getNumEdgeNode(*it);
}

inline void count_edge_3()
{
    size_t cnt = num_refresh;
    while (cnt > 0)
    {
        // cout << g_hg.getEdge(cnt - 1).size() <<  "\t";
        if (g_hg.getEdge(cnt - 1).size() >= 2)
        {
            g_edgeList3.push_back(g_hg.getEdge(cnt - 1));
        }
        cnt--;
    }
    // cout << endl;
}

#endif