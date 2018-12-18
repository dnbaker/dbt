#include "pfparse.h"
#include <fstream>
#include <iostream>
#include <omp.h>
#include "klib/ketopt.h"

using namespace dbt;

void usage() {
    std::fprintf(stderr, "Usage: ctest <genome1> <genome2> ...\nPerforms transform on genome.\n");
    std::exit(EXIT_FAILURE);
}

// static int ketopt(ketopt_t *s, int argc, char *argv[], int permute, const char *ostr, const ko_longopt_t *longopts)

std::vector<std::string> load_from_file(const char *path) {
    std::vector<std::string> ret;
    std::ifstream ifs(path);
    std::string line;
    while(std::getline(ifs, line))
        ret.emplace_back(line);
    return ret;
}

int main(int argc, char *argv[]) {
    omp_set_num_threads(1);
    size_t wsz = 10;
    unsigned nthreads = 1;
    ketopt_t o = KETOPT_INIT;
    std::vector<std::string> paths;
    for(int c; (c = ketopt(&o, argc, argv, true /* permute */, "w:p:F:h?", nullptr)) >= 0;) {
        switch(c) {
            case 'p': nthreads = std::atoi(o.arg); break;
            case 'F': paths = std::move(load_from_file(o.arg)); break;
            case 'w': wsz = std::atoi(o.arg); break;
        }
    }
    dbt::HashPasser<std::FILE *> hp(argv[o.ind], wsz, nthreads, true);
    hp.run();
}
