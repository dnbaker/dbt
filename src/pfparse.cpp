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



std::string sanitize(std::string in) {                                                                      
    if(auto i = in.rfind('/'); i != std::string::npos) in = std::string(i + 1 + in.begin(), in.end()); 
    return in;                                                                                         
}    

int main(int argc, char *argv[]) {
    if(std::find_if(argv, argv + argc, [](auto x) {return std::strcmp(x, "--help") == 0;}) != argv + argc)
        usage();
    size_t wsz = 10;
    unsigned nthreads = 1;
    ketopt_t o = KETOPT_INIT;
    std::vector<std::string> paths;
    std::string final_prefix = "final_prefix";
    for(int c; (c = ketopt(&o, argc, argv, true /* permute */, "w:p:F:h?", nullptr)) >= 0;) {
        switch(c) {
            case 'O': final_prefix = o.arg; break;
            case 'p': nthreads = std::atoi(o.arg); break;
            case 'F': paths = std::move(load_from_file(o.arg)); break;
            case 'w': wsz = std::atoi(o.arg); break;
            case 'h': case '?': usage();
        }
    }
    if(nthreads > 1000) throw "a party";
    omp_set_num_threads(nthreads);
    if(argc == o.ind && paths.empty()) usage(); // No paths
    if(paths.empty())
        paths = std::vector<std::string>(argv + o.ind, argv + argc);
    #pragma omp parallel for schedule(static,1)
    for(size_t i = 0; i < paths.size(); ++i) {
        dbt::HashPasser<gzFile> hp(paths[i].data(), wsz, 1, true);
        hp.run(sanitize(paths[i]).data());
    }
    #pragma omp parallel for schedule(static,1)
    for(size_t i = 0; i < paths.size(); ++i) {
        paths[i] = sanitize(paths[i]);
    }
    omp_set_num_threads(1); // Return to OS.
    massive_merge(paths, final_prefix, nthreads);
}
