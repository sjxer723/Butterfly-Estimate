#include "btf-count.h"

Graph::Graph(string path, config_t c) {
    measure_query = c.measure_query;
    measure_time = c.measure_time;

    printf("Loading graph...\n");
    this->is_original = is_original;
    this->is_bipartite = 0;
    FILE *fin = nullptr;
    std::string binfile_name =
        "graph-sort" + ((c.density != 10) ? ("-0" + to_string(c.density)) : "") + ".bin";

    path = append_dir(path, binfile_name);
    fin = fopen(path.c_str(), "rb");
    if (!fin) {
        fprintf(stderr, "Fail to open file %s!\n", path.c_str());
        exit(0);
    }
    size_t items_read = 0;

    items_read = fread(&lsize, sizeof(int), 1, fin);
    items_read = fread(&n, sizeof(int), 1, fin);
    items_read = fread(&m, sizeof(ll), 1, fin);
    deg = new int[n];
    deg_cnt = new int[n];
    dat = new int[m];
    con = new int *[n];
    nq = new bool[n];
    mq = new bool[m];
    logn = calc_log(n);
    loglogn = calc_log(logn);

    items_read = fread(deg, sizeof(int), this->n, fin);
    if (items_read != n) {
        printf("Failed to read node data.\n");
        fclose(fin);
        return;
    }

    items_read = fread(dat, sizeof(int), this->m, fin);
    fclose(fin);
    memset(mq, 0, sizeof(bool) * m);
    memset(nq, 0, sizeof(bool) * n);
    m /= 2; // undirect
    theta = (int)(sqrt(6.0 * m));
    if ((ll)(theta)*theta < m)
        ++theta;
    m_bar = m; // temprarily use exact m,w to replace m_bar and w_bar, need to
               // be improved
    // wedge_cnt();

    ll p = 0;
    for (int i = 0; i < n; ++i) {
        con[i] = dat + p;
        p += deg[i];
        deg_cnt[i] = p;
    }

    printf("%s:\n", path.c_str());
    printf("deg_cnt[n-1] = %d\n", deg_cnt[n - 1]);
    printf("deg[lsize-1] = %d\n", deg_cnt[lsize - 1]);
    printf("n1 = %d\n", lsize);
    printf("n=%d,m=%lld\n", n, m);
    to_bip = c.uni_to_bip;
    alpha = 1;
    query_cnt = -1;
    BTF_cnt = -1;
    BTF_time = -1;
}

Graph::~Graph() {
    delete[] deg;
    delete[] dat;
    delete[] con;
    // if (cnt) delete[] cnt; if (cnt_dat) delete[] cnt_dat;
}

void Graph::save_expr_results(string res_file_name) {
    FILE *res_file = fopen(res_file_name.c_str(), "w");

    printf("BTF: %lld,\tTime: %.3f,\tQuery: %lld\n", BTF_cnt, BTF_time, query_cnt);
    fprintf(res_file, "n:%d\nm:%lld\nBTF:%lld\nTime: %.3f\nQuery: %lld\n", n, m, BTF_cnt, BTF_time,
            query_cnt);
}

void Graph::rec_count_arb() {
    double dlt = 0;
    long long cnt_rec = 0, cnt_tmp = 0, cnt_bistar = 0;
#pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads(), u, v, w;

        long long local_cnt = 0, local_cnt_tmp = 0, local_cnt_bistar = 0;
        double local_dlt = 0;
        int *last_use = new int[n], *last_cnt = new int[n], pre = 0;
        memset(last_use, -1, sizeof(int) * n);
        int idx = 0;
        for (u = pid; u < n; u += np) {
            for (int j = 0; j < idx; j++) {
                last_cnt[last_use[j]] = 0;
            }
            idx = 0;
            for (int i = 0; i < query_deg(u); ++i) {
                v = query_nbr(u, i);
                pre = -1;
                for (int j = 0; j < query_deg(v); ++j) {
                    w = query_nbr(v, j);
                    ++local_cnt_tmp;
                    if (w >= u || w >= v)
                        break;
                    ++local_cnt_bistar;
                    if (pre == -1)
                        pre = w;
                    else {
                        local_dlt += abs(w - pre);
                        pre = w;
                    }
                    local_cnt += last_cnt[w];
                    ++last_cnt[w];
                    if (last_cnt[w] == 1)
                        last_use[idx++] = w;
                }
            }
        }
#pragma omp critical
        {
            cnt_rec += local_cnt;
            cnt_tmp += local_cnt_tmp;
            cnt_bistar += local_cnt_bistar;
            dlt += local_dlt;
        }
    }
    BTF_cnt = cnt_rec;
}

void Graph::btf_count_ESpar(double p) {
    m = 0;
    vector<int> *new_con = new vector<int>[MAXN];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < deg[i]; ++j) {
            if (con[i][j] < i)
                continue;
            if (p * RAND_MAX > rand()) {
                new_con[i].emplace_back(con[i][j]);
                new_con[con[i][j]].emplace_back(i);
                m += 2;
            }
        }
    }
    delete dat;
    dat = new int[m];
    int pnt = 0;
    for (int i = 0; i < n; ++i) {
        con[i] = dat + pnt;
        deg[i] = new_con[i].size();
        for (int j = 0; j < deg[i]; ++j)
            dat[pnt++] = new_con[i][j];
    }
    m /= 2;
    this->rec_count_arb();
    if (BTF_cnt == 0)
        return;
    BTF_cnt = BTF_cnt / pow(p, 4);
}

/**
 * sample_wedge(e, seed): Sample a wedge from an input edge
 * @param e: the input edge (u, v)
 * @param seed: the random number seed
 */
Graph::wedge Graph::sample_wedge(edge e, int *seed, size_t deg_u, size_t deg_v) {
    // uniformly sample an additional node from the neighbors of `u` and `v`
    int idx = randnum(deg_u + deg_v, seed);
    int w = (idx < deg_u ? query_nbr(e.u, idx) : query_nbr(e.v, idx - deg_u));
    if (w != e.u && w != e.v)
        return (idx < deg_u ? wedge(e.u, e.v, w) : wedge(e.v, e.u, w));
    else
        return wedge(-1, -1, -1);
}

bool Graph::is_butterfly(wedge w, int z) {
    if (z == w.w || z == w.u || z == w.v)
        return false;
    return query_pair(w.u, z) && query_pair(w.v, z);
}

ll Graph::btf_count_est(int s1, int s2) {
    int thread_n, id, org_seed = rand();
    int seed = org_seed;
    ll sum_deg = 0;
    vector<ll> cnt;
    vector<edge> S;
    vector<pair<size_t, size_t>> S_degree;

    // Sample s1 edges uniformly
    for (int i = 0; i < s1; i++) {
        size_t edge_idx = randnum(2 * m, &seed);
        edge e = query_edge(edge_idx);
        size_t deg_u = query_deg(e.u), deg_v = query_deg(e.v);
        sum_deg += deg_u + deg_v - 2;
        cnt.emplace_back(sum_deg);
        S.emplace_back(e); // store the sampled edges
        S_degree.emplace_back(
            std::make_pair(deg_u, deg_v)); // store the degrees of these edge pairs
    }

    double cnt_x = 0;
#pragma omp parallel private(id, seed) shared(thread_n)
    {
        thread_n = omp_get_num_threads();
        id = omp_get_thread_num();
        seed = org_seed + id * 256;
        ll r_sum = 0;
#pragma omp for reduction(+ : cnt_x)
        for (int i = 0; i < s2; ++i) {
            ll tmp = randnum(sum_deg, &seed);
            int idx = upper_bound(cnt.begin(), cnt.end(), tmp) - cnt.begin();
            edge e = S[idx];
            std::pair<size_t, size_t> deg_uv = S_degree[idx];
            wedge wdg = sample_wedge(e, &seed, deg_uv.first,
                                     deg_uv.second); // return wdg(-1,-1,-1) when sample fails.
            if (wdg.u == -1)
                continue;
            int y = wdg.u, z = wdg.v;
            int deg_z = query_deg(z);
            int deg_y = (y == e.u) ? deg_uv.first : deg_uv.second;
            int deg_v = deg_z; // v is the sampled node
            int deg_y_record = deg_y;
            if (smaller(z, deg_z, y, deg_y)) {
                swap(y, z);
                swap(deg_y, deg_z);
            }
            double root_m = sqrt(m_bar * 1.0);
            int r = ceil(deg_y / root_m);
            if (deg_y <= root_m) {
                tmp = randnum(m_bar, &seed);
                if (tmp > deg_y * root_m)
                    r = 0;
            }
            double cnt_y = 0;
            if (r > 0) {
                for (int j = 0; j < r; ++j) {
                    int x = query_nbr(y, randnum(deg_y, &seed)), deg_x = query_deg(x);
                    if (is_butterfly(wdg, x) && smaller(wdg.v, deg_v, x, deg_x)) {
                        cnt_y += max(root_m, (double)(deg_y)) * 0.25;
                    }
                }
                cnt_y /= r;
            }
            r_sum += (ll)r;
            cnt_x += cnt_y;
        }
    }
    ll num_of_btf = (cnt_x * m * sum_deg / S.size() / s2);

    return num_of_btf;
}

void Graph::estimate(double s1_coeff, double s2_coeff, int r_base, int r_exp, int r_round) {
    ll res = 0;
    int s_cnt = 0;
    int s1 = sqrt(m) * s1_coeff, s2 = sqrt(m) * s2_coeff;
    int s_tmp_base = r_base, s_tmp_exp = r_exp, s_tmp = r_base; // run base * exp^i rounds
    float error_rate = 0;
    ll ground_truth = 0;

    for (int i = 0; i < r_round + 1; i++) {
        int to_run_rounds = (i > 0) ? (s_tmp_exp - 1) * s_tmp : s_tmp;
        res = (i > 0) ? res / (s_tmp_exp) : res;
        s_tmp = (i > 0) ? s_tmp * (s_tmp_exp) : s_tmp;
// Get the average number of `s_tmp` estimation
#pragma omp parallel for reduction(+ : res)
        for (int k = 0; k < to_run_rounds; k++) {
            ll tmp = btf_count_est(s1, s2);
            res += tmp / s_tmp;
        }
        if (ground_truth > 0) {
            error_rate = abs((res - ground_truth) / (ground_truth * 1.0));
            printf("%d: %f\n", s_tmp, error_rate);
        }
    }
    BTF_cnt = res;
}

void Graph::weighted_one_side_sampling(int r_base, int r_exp, int r_round) {
    ll deg_count = 0;
    ll btf_num = 0;
    vector<int> deg_pos_vertex_vec;
    vector<ll> deg_cnt_vec;
    int round_num_base = r_base, round_num = round_num_base;
    chrono::duration<double> elapsed_seconds;

    // If the bipartite graph is constructed from a unipartite graph,
    // L is the set of vertices with even indeices
    for (int v = 0; v < to_bip ? n : lsize; v += ((to_bip) ? 2 : 1)) {
        int v_deg = query_deg(v);
        if (v_deg) {
            deg_count += v_deg;
            deg_cnt_vec.emplace_back(deg_count);
            deg_pos_vertex_vec.emplace_back(v);
        }
    }
    printf("Sum of deg(v):%lld\n", deg_count);
    printf("m:%lld\n", m);

    // for (int i=0; elapsed_seconds.count() <= 64; i++) {
    for (int i = 0; i < r_round; i++) {
        int to_run_rounds = (i > 0) ? (r_exp - 1) * round_num : round_num;
        btf_num = (i > 0) ? btf_num / r_exp : btf_num;
        round_num = (i > 0) ? round_num * r_exp : round_num;

#pragma omp parallel for reduction(+ : btf_num)
        for (int r = 0; r < to_run_rounds; r++) {
            int seed1 = rand(), seed2 = rand();
            int u = 0, v = 0;
            int xi_uv = 0;
            ll tmp1 = randnum(deg_count, &seed1), tmp2 = randnum(deg_count, &seed2);

            u = deg_pos_vertex_vec[upper_bound(deg_cnt_vec.begin(), deg_cnt_vec.end(), tmp1) -
                                   deg_cnt_vec.begin()];
            v = deg_pos_vertex_vec[upper_bound(deg_cnt_vec.begin(), deg_cnt_vec.end(), tmp2) -
                                   deg_cnt_vec.begin()];
            if (u == v) {
                continue;
            }
            int u_deg = query_deg(u), v_deg = query_deg(v);
            if (u_deg > 0 && v_deg > 0) {
                if (u_deg > v_deg) {
                    swap(u, v);
                }
                for (int w_idx = 0, u_deg = query_deg(u); w_idx < u_deg; w_idx++) {
                    int w = query_nbr(u, w_idx);
                    if (query_pair(v, w)) {
                        xi_uv++;
                    }
                }
                ll xi_uv_choose2 = (ll)xi_uv * ((ll)xi_uv - 1) / 2;
                ll tmp_btf = (m / (2 * u_deg) * (m / v_deg)) * xi_uv_choose2;
                btf_num += tmp_btf / round_num;
            }
        }
    }
    BTF_cnt = btf_num;
}

void btf_count(config_t c, int method) {
    if (c.n_threads > 0)
        omp_set_num_threads(c.n_threads);
    Graph *g = new Graph(c.path, c);
    string output_filename;
    ;

    auto start = chrono::system_clock::now();
    if (method == 0) {
        // Run exact butterfly counting algorithm
        output_filename = append_dir(c.path, "Exact_res.txt");
        g->rec_count_arb();
    } else if (method == 1) {
        // Run ESpar butterfly counting algorithm
        output_filename = append_dir(c.path, "Espar_res.txt");
        g->btf_count_ESpar(0.2);
    } else if (method == 2) {
        // Run TLS(Two Leveled Sampling) butterfly counting algorithm
        output_filename = append_dir(c.path, "Tls_res.txt");
        g->estimate(c.s1, c.s2, c.r_base, c.r_exp, c.r_round);
    } else if (method == 3) {
        // Run WPS(Weighted one-sided Pair Sampling) butterfly counting algorithm
        output_filename = append_dir(c.path, "Wps_res.txt");
        g->weighted_one_side_sampling(c.r_base, c.r_exp, c.r_round);
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> estimate_time = end - start;
    g->BTF_time = estimate_time.count();

    // Save the experimental results
    g->save_expr_results(output_filename);
    delete g;
}