#include "pfparse.h"

using namespace dbt;

void usage() {
    std::fprintf(stderr, "Usage: ctest <genome1> <genome2> ...\nPerforms transform on genome.\n");
    std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    size_t wsz = 30;
    std::vector<dbt::HashPass<>> vecs;
    dbt::khmap main_map;
    std::vector<std::string> sans;
    if(argc == 1) usage();
    for(auto av = argv + 1; *av; ++av) {
        std::string v = *av, v2 = v;
        std::string san(util::sanitize_name(v));
        //if(auto l = v2.find('.'); l != std::string::npos && l) v2.resize(l);
        std::fprintf(stderr, "About to create thing at window size %zu with input path %s (v2 data: %s)\n", wsz, *av, v2.data());
        vecs.emplace_back(wsz, 1, 0, san.data());
        auto &h = vecs.back();
        h.open_ifp(v.data());
        h.make_map();
        h.fill();
        //std::fprintf(stderr, h.str().data());
        h.map_->assert_nonnull();
        sans.emplace_back(san);
    }
    std::fprintf(stderr, "Now merging stuff together. memory usage: %zu\n", vecs[0].memory_usage());
}
