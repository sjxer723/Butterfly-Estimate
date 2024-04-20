#ifndef BTFCOUNT_H_
#define BTFCOUNT_H_

#include "utility.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

#define MAXN 2000000000

// Butterfly counting algorithms
void btf_count(config_t c, int method);

class Graph {
  private:
    ll theta; // sqrt(6m)
    ll w_bar; // estimation of number of wedges
    ll b_bar, b_tld;
    char is_bipartite; // whether the graph is a bipartitie
    bool *mq, *nq;
    double alpha;          // tuning the size of s1 & s2
    int lsize;             // size of left vertex set
    int n, logn, loglogn;  // number of vertices
    int *dat, **con, *deg; // data, edge, degree
    int *deg_cnt;
    bool *dat_mark, *deg_mark;
    // char *node_type;
    bool is_original;
    bool measure_query;
    bool measure_time;

    struct edge {
        long long u, v;
        edge() {}
        edge(int u, int v) : u(u), v(v) {}
    };

    /**
     * Wedge is a structure used to represent a pair of edges.
     * i.e., a wedge (w, u, v) represent edge (w, u) and (w, v), which
     * share the same vertex w.
     */
    struct wedge {
        long long w, u, v;
        wedge() {}
        wedge(int w, int u, int v) : w(w), u(u), v(v) {}
    };

  public:
    ll m;     // number of edges
    ll m_bar; // estimation of m
    ll BTF_cnt;
    double BTF_time; // runtime of alg
    ll query_cnt;    // counter for queries
    bool to_bip;     // indicates whether the graph is constructed from a unipartite graph
    void wedge_cnt();

    // Allowed query
    inline int query_deg(int v);
    inline bool query_pair(int v, int u);
    inline edge query_edge(size_t idx);
    inline int query_nbr(int v, int i);

    struct wedge sample_wedge(edge e, int *seed, size_t deg_u, size_t deg_v);
    bool is_butterfly(wedge w, int z);
    ll btf_count_est(int s1, int s2);

    // Butterfly counting methods
    void rec_count_arb();                                                                // exact
    void btf_count_ESpar(double p);                                                      // espar
    void estimate(double s1_coeff, double s2_coeff, int r_base, int r_exp, int r_round); // tls
    void weighted_one_side_sampling(int r_base, int r_exp,
                                    int r_round); // weighted one-sided pair sampling
    void save_expr_results(string res_file_name);

  public:
    Graph(string path, config_t c);
    ~Graph();

}; // Class Graph

inline int Graph::query_deg(int v) {
    if (measure_query) {
#pragma omp atomic
        ++query_cnt;
    }

    return deg[v];
}

inline bool Graph::query_pair(int v, int u) {
    if (measure_query) {
#pragma omp atomic
        ++query_cnt;
    }
    int *tmp = lower_bound(con[v], con[v] + deg[v], u);
    char query_res = (tmp != con[v] + deg[v]) && *tmp == u;

    return query_res;
}

inline int Graph::query_nbr(int v, int i) {
    if (measure_query) {
#pragma omp atomic
        ++query_cnt;
    }

    return con[v][i];
}

inline Graph::edge Graph::query_edge(size_t idx) {
    if (measure_query) {
#pragma omp atomic
        ++query_cnt;
    }
    size_t left = 0;
    size_t right = n - 1;
    size_t u = 0;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (idx <= deg_cnt[mid]) {
            if (mid == 0 || idx > deg_cnt[mid - 1]) {
                u = mid;
                break;
            } else {
                right = mid - 1;
            }
        } else {
            left = mid + 1;
        }
    }
    assert(u >= 0 && idx <= deg_cnt[u]);
    if (u == 0) {
        return edge(u, con[u][idx]);
    } else {
        return edge(u, con[u][idx - deg_cnt[u - 1]]);
    }
}

inline void Graph::wedge_cnt() {
    w_bar = 0;
    for (int i = 0; i < n; ++i)
        w_bar += (ll)deg[i] * (ll)(deg[i] - 1) / 2;
}

#endif