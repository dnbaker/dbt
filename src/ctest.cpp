#include "pfparse.h"

int main() {
    dbt::HashPass<> hp(1337);
    dbt::HashPass<gzFile> h2p(1337);
    std::vector<dbt::HashPass<>> vecs;
    for(size_t i = 0; i < 20; vecs.push_back(1337));
    dbt::merge_hashpasses("sdfds", vecs, nullptr);
}
