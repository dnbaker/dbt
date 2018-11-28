#ifndef UTIL_H____
#define UTIL_H____
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif
#include <cstdint>
#include <cstdlib>
#include <cstring>

#ifndef FOREVER
#  define FOREVER for(;;)
#endif
namespace dbt {
namespace util {

} // util
} // dbt

#endif // #ifndef UTIL_H____
