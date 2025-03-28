#ifndef UTILITY_H
#define UTILITY_H

#include <algorithm>
#include <climits>
#include <cstring>
#include <string>
#include <vector>

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
#ifndef MAXN
#define MAXN 2000000000
#endif
#ifndef MAXN1
#define MAXN1 1000000000
#endif
#define INC_FACTOR 1.2

typedef struct Config {
    string path;
    string method;
    int density = 10;
    int n_threads = 1;
    bool measure_query = false;
    bool measure_time = true;
    bool is_bipartite_format = true;
    bool uni_to_bip = false;
    int r_base = 50;
    int r_exp = 2;
    int r_round = 0;
    float s1 = 1;
    float s2 = 1;
} config_t;

class IntNode {
  public:
    int val, id;
    bool operator<(const IntNode &v) const {
        if (val == v.val)
            return id < v.id;
        else
            return val < v.val;
    };
};

inline string append_dir(string path, string filename) {
    string new_dir;
    if (path.back() == '/') {
        new_dir = path + filename;
    } else {
        new_dir = path + "/" + filename;
    }

    return new_dir;
}

inline int random_value(int *seed) {
#ifndef MONTECARLO
    return rand();
#endif
    *seed = ((3881 * (*seed) + 509) % 32768);
    return *seed;
}

inline ll rand64(int *seed) {
    ll ret = ((((ll(random_value(seed))) << 60) + ((ll(random_value(seed))) << 45) +
               ((ll(random_value(seed))) << 30) + ((ll(random_value(seed))) << 15) +
               ll(random_value(seed))) &
              (LLONG_MAX >> 1));
    return ret;
}

inline int rand32() {
    return (((ll(rand())) << 30) + ((ll(rand())) << 15) + ll(rand())) & (INT_MAX >> 1);
}

inline ll randnum(ll range, int *seed) {
    ll rand64_val = rand64(seed);
    ll ret = rand64_val % range;
    return ret;
}

inline bool smaller(int a, int a_deg, int b, int b_deg) {
    return a_deg < b_deg || (a_deg == b_deg && a < b);
}

inline int calc_log(int num) {
    int cnt = 0;
    while (num > 0) {
        ++cnt;
        num /= 2;
    }
    return cnt;
}

// bern(): return a Bernoulli variable with 1's probability p
inline int bern(int p_nmr, int p_dnr, int *seed) {
    int tmp = randnum(p_dnr, seed);
    return (tmp <= p_nmr ? 1 : 0);
}

inline int mygcd(const int &aa, const int &bb) {
    int a = aa, b = bb;
    int r;
    while (b > 0) {
        r = a % b;
        a = b;
        b = r;
    }
    return a;
}

inline void reduce(int &nmr, int &dnr) {
    int tmp = mygcd(nmr, dnr);
    if (tmp > 1) {
        nmr /= tmp;
        dnr /= tmp;
    }
}

inline int get_num_cnt(string path) {
    FILE *fin = fopen((path + "graph.txt").c_str(), "r");
    char line[1000], st[64];
    int cnt = 0, min_cnt = 100;

    while (fgets(line, 1000, fin) && cnt < 10) {
        if (!isdigit(line[0]))
            continue;
        vector<char *> v_num;
        int len = (int)strlen(line);
        for (int i = 0; i < len; ++i)
            if (!isdigit(line[i]) && line[i] != '.')
                line[i] = '\0';
            else if (i == 0 || !line[i - 1])
                v_num.push_back(line + i);
        if ((int)v_num.size() < 2)
            continue;
        min_cnt = min(min_cnt, (int)v_num.size());
        ++cnt;
    }
    fclose(fin);
    return min_cnt;
}

inline bool get_edge(char *line, int &a, int &b, int num_cnt) {
    if (!isdigit(line[0]))
        return false; // for other datasets
    vector<char *> v_num;
    int len = (int)strlen(line);
    for (int i = 0; i < len; ++i)
        if (!isdigit(line[i]) && line[i] != '.')
            line[i] = '\0';
        else if (i == 0 || !line[i - 1])
            v_num.push_back(line + i);
    if ((int)v_num.size() != num_cnt)
        return false;
    sscanf(v_num[0], "%d", &a);
    sscanf(v_num[1], "%d", &b);
    return true;
}

inline void get_order(vector<int> *con, int n, int *o) {
    IntNode *f = new IntNode[n];
    for (int i = 0; i < n; ++i)
        f[i].id = i, f[i].val = (int)con[i].size();
    sort(f, f + n);
    for (int i = 0; i < n; ++i)
        o[i] = f[n - 1 - i].id;
    delete[] f;
}

#endif