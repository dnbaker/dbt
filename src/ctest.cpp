#include "pfparse.h"

using namespace dbt;

void usage() {
    std::fprintf(stderr, "Usage: ctest <genome1> <genome2> ...\nPerforms transform on genome.\n");
    std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    size_t wsz = 200;
    std::vector<dbt::HashPass<>> vecs;
    dbt::khmap main_map;
    if(argc == 1) usage();
    for(auto av = argv + 1; *av; ++av) {
        std::string v = *av, v2 = v;
        if(auto l = v2.find('.'); l != std::string::npos && l) v2.resize(l);
        vecs.emplace_back(wsz, 1, 0, v2.data());
        auto &h = vecs.back();
        h.open_ifp(v.data());
    }
    dbt::merge_hashpasses("sdfds", vecs, &main_map);
}
