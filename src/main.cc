#include "btf-count.h"
#include "utility.h"
#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <unistd.h>

using namespace std;

config_t c;

void txt_to_bin(config_t c) {
    FILE *fin;
    fin = fopen((c.path + "graph.txt").c_str(), "r");
    char line[1024];
    int n = 0, a, b, n1 = 0, n2 = 0, num_cnt = get_num_cnt(c.path);
    vector<int> *con = new vector<int>[MAXN];
    long long cnt = 0, m = 0;

    printf("Loading text, num_cnt=%d...\n", num_cnt);
    while (fgets(line, 1024, fin)) {
        if (!get_edge(line, a, b, num_cnt))
            continue;
        if (a < 0 || b < 0)
            continue;
        if (c.uni_to_bip && ((a & 1) == (b & 1))) {
            continue;
        }
        if (rand() % 10 >= c.density) {
            continue;
        }
        // printf("a = %d, b= %d\n", a, b);
        n1 = max(n1, a + 1);
        n2 = max(n2, b + 1);
        if (c.is_bipartite_format)
            b += MAXN1;
        else if (a == b)
            continue;
        if (con[a].capacity() == con[a].size())
            con[a].reserve(con[a].size() * INC_FACTOR);
        if (con[b].capacity() == con[b].size())
            con[b].reserve(con[b].size() * INC_FACTOR);
        con[a].push_back(b);
        con[b].push_back(a);
        if ((++cnt) % (long long)10000000 == 0)
            printf("%lldM lines finished\n", cnt / 1000000);
    }
    fclose(fin);
    n = c.is_bipartite_format ? n1 + n2 : max(n1, n2);

    printf("n1=%d,n2=%d\n", n1, n2);

    if (c.is_bipartite_format) {
        printf("Re-numbering bipartite...\n");
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < (int)con[i].size(); ++j)
                con[i][j] = con[i][j] - MAXN1 + n1;
        for (int i = 0; i < n2; ++i) {
            con[n1 + i].clear();
            con[n1 + i].reserve(con[MAXN1 + i].size());
            for (int j = 0; j < (int)con[MAXN1 + i].size(); ++j)
                con[n1 + i].push_back(con[MAXN1 + i][j]);
            vector<int>().swap(con[MAXN1 + i]);
        }
    }

    int maxd = 0;
    for (int i = 0; i < n; ++i)
        if (con[i].size() > 0) {
            sort(con[i].begin(), con[i].end());
            int p = 1;
            for (int j = 1; j < (int)con[i].size(); ++j)
                if (con[i][j - 1] != con[i][j])
                    con[i][p++] = con[i][j];
            con[i].resize(p);
            m += p;
            maxd = max(maxd, p);
        }

    printf("Saving binary...\n");
    int *deg = new int[n];
    for (int i = 0; i < n; ++i)
        deg[i] = (int)con[i].size();

    FILE *fout = fopen(
        (c.path + "graph-sort" + ((c.density != 10) ? ("-0" + to_string(c.density)) : "") + ".bin")
            .c_str(),
        "wb");
    fwrite(&n1, sizeof(int), 1, fout);
    fwrite(&n, sizeof(int), 1, fout);
    fwrite(&m, sizeof(long long), 1, fout);
    fwrite(deg, sizeof(int), n, fout);

    int *nbr = new int[maxd];
    for (int i = 0; i < n; ++i) {
        int d = con[i].size();
        for (int j = 0; j < d; ++j)
            nbr[j] = con[i][j];
        sort(nbr, nbr + d);
        fwrite(nbr, sizeof(int), d, fout);
    }

    fclose(fout);
    printf("Created binary file, n = %d, m = %lld\n", n, m);
    delete[] con;
    delete[] deg;
}

void get_config(int argc, char **argv) {
    int opt;
    int option_index = 0;

    static struct option long_options[] = {
        {"rbase", required_argument, NULL, 'b'},  {"rexp", required_argument, NULL, 'e'},
        {"rround", required_argument, NULL, 'r'}, {"s1", required_argument, NULL, 'o'},
        {"s2", required_argument, NULL, 'v'},     {"uni-to-bip", no_argument, NULL, 'u'},
        {"is-uni", no_argument, NULL, 'i'},       {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "p:m:n:d:qt", long_options, &option_index)) != -1) {
        switch (opt) {
        case 'b':
            c.r_base = atoi(optarg);
            break;
        case 'e':
            c.r_exp = atoi(optarg);
            break;
        case 'r':
            c.r_round = atoi(optarg);
            break;
        case 'o':
            c.s1 = atof(optarg);
            break;
        case 'v':
            c.s2 = atof(optarg);
            break;
        case 'i':
            c.is_bipartite_format = false;
            break;
        case 'u':
            c.uni_to_bip = true;
            break;
        case 'p':
            c.path = optarg;
            printf("input graph path: %s\n", c.path.c_str());
            break;
        case 'm':
            c.method = optarg;
            printf("method: %s\n", c.method.c_str());
            break;
        case 'n':
            c.n_threads = atoi(optarg);
            break;
        case 'd':
            c.density = atoi(optarg);
            break;
        case 'q':
            c.measure_query = true;
            break;
        case 't':
            c.measure_time = true;
            break;
        case '?': // unknown option...
            cerr << "Unknown option: '" << char(optopt) << "'!" << endl;
            break;
        }
    }
}

int main(int argc, char **argv) {
    get_config(argc, argv);

    srand((int)time(0));
    if (strcmp(c.method.c_str(), "exact") == 0) {
        btf_count(c, 0);
    } else if (strcmp(c.method.c_str(), "espar") == 0) {
        btf_count(c, 1);
    } else if (strcmp(c.method.c_str(), "tls") == 0) {
        btf_count(c, 2);
    } else if (strcmp(c.method.c_str(), "wps") == 0) {
        btf_count(c, 3);
    } else if (strcmp(c.method.c_str(), "txt2bin") == 0) {
        txt_to_bin(c);
    }

    return 0;
}