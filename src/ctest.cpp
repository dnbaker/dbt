#include "pfparse.h"

using namespace dbt;

void usage() {
    std::fprintf(stderr, "Usage: ctest <genome1> <genome2> ...\nPerforms transform on genome.\n");
    std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    size_t wsz = 10;
    dbt::HashPasser<> hp(argv[1], 10, 1, true);
    hp.run();
}
