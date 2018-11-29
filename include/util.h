#ifndef UTIL_H____
#define UTIL_H____
#include <cstdint>
#include <cstring>
#include "fpwrap.h"

#ifndef FOREVER
#  define FOREVER for(;;)
#endif
namespace dbt {
namespace util {
using namespace ::fp;

char *strdup(const char *s) {
    size_t l = std::strlen(s) + 1;
    char *ret = static_cast<char *>(std::malloc(l));
    if(!ret) throw std::bad_alloc();
    std::memcpy(ret, s, l);
    return ret;
}


} // util
} // dbt

#endif // #ifndef UTIL_H____
