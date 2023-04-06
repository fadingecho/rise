#ifndef FITNESS_FUNCTIONS_H
#define FITNESS_FUNCTIONS_H

#include <vector>
#include <iostream>
#include <memory>
#include <set>

#include "global.hpp"

extern Graph g_g;
extern HyperGraph g_hg;

class IMP
{
public:
    explicit IMP(unsigned int n, int t, bool mo)
    {
        n_ = n;
        mo_ = mo;
    }

    vector<double> operator()(const vector<int> &x)
    {
        extern vector<int> g_nodeDegree, g_nodeDegree_cpy;
        extern vector<bool> g_edgeMark;

        // All objective functions are to be maximized using the GAs.
        // Note: All of the functions are modified for maximization where needed, and return vector<double>.
        int inf = 0;

        for (vector<int>::const_iterator it = g_hg.getNodeMapVec().cbegin(); it != g_hg.getNodeMapVec().cend(); it++)
            g_nodeDegree_cpy[*it - 1] = g_nodeDegree[*it - 1];

        size_t numEdge = g_hg.getNumEdge();

        // check if an edge is removed
        g_edgeMark.assign(numEdge, false);

        // calculate inf
        // note that x.size == n, node idx starts from 1
        vector<int>::const_iterator n_it, e_it, n_it_end, e_it_end;
        for (vector<int>::const_iterator it = x.cbegin(); it != x.cend(); it++)
        {
            if (!g_hg.testNodeMap(*it + 1))
                continue;

            inf += g_nodeDegree_cpy[*it];
            e_it_end = g_hg.getNode(*it + 1).cend();
            for (e_it = g_hg.getNode(*it + 1).cbegin(); e_it != e_it_end; ++e_it)
            {
                int e_i = *e_it;
                if (g_edgeMark[e_i])
                {
                    continue;
                }
                
                n_it_end = g_hg.getEdge(e_i).cend();
                for (n_it = g_hg.getEdge(e_i).cbegin(); n_it != n_it_end; ++n_it)
                {
                    g_nodeDegree_cpy[*n_it - 1]--;
                }
                g_edgeMark[e_i] = true;
            }
        }

        return {static_cast<double>(inf) / static_cast<double>(g_hg.getNumEdge()) * static_cast<double>(n_), 1.0};
    }

    inline bool getMo()
    {
        return mo_;
    }

private:
    unsigned int n_;
    bool mo_;
};

#endif // !FITNESS_FUNCTIONS_H