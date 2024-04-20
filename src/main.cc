#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <unistd.h>

using namespace std;

struct Config {
    string path;
    string method;
    bool flagA = false;
    bool flagB = false;
};

struct Config get_config(int argc, char **argv) {
    int opt;
    struct Config c;

    while ((opt = getopt(argc, argv, "p:m:")) != -1) {
        switch (opt) {
        case 'p':
            c.path = optarg;
            printf("input graph path: %s\n", c.path.c_str());
            break;
        case 'm':
            c.method = string(optarg);
            printf("method: %s\n", c.path.c_str());
            break;
        case '?': // unknown option...
            cerr << "Unknown option: '" << char(optopt) << "'!" << endl;
            break;
        }
    }

    return c;
}

int main(int argc, char **argv) {
    get_config(argc, argv);

    return 0;
}