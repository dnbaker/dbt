#ifndef UTIL_H____
#define UTIL_H____
#include <cstdint>
#include <cstring>
#include "fpwrap/fpwrap.h"
#include <cinttypes>

#ifndef FOREVER
#  define FOREVER for(;;)
#endif

#if !defined(unlikely)
#  if __GNUC__ || __clang__
#    define unlikely(x) __builtin_expect((x), 0)
#  else
#    define unlikely(x) (x)
#  endif
#endif
#if !defined(likely)
#  if __GNUC__ || __clang__
#    define likely(x) __builtin_expect((x), 1)
#  else
#    define likely(x) (x)
#  endif
#endif

namespace dbt {
namespace util {
using namespace ::fp;
enum SpecialChars: std::uint8_t {Dollar = 0x01, EndOfWord = 0x00};

static inline char *strdup(const char *s) {
    size_t l = std::strlen(s) + 1;
    char *ret = static_cast<char *>(std::malloc(l));
    if(!ret) throw std::bad_alloc();
    std::memcpy(ret, s, l);
    return ret;
}

static std::string sanitize_name(std::string n) {
    auto s = n.find_last_of('/');
    s = s == std::string::npos ? 0: s;
    auto p = n.find('.', s);
    p = p == std::string::npos ? n.size(): p;
    return n.substr(s, p - s);
}


} // util
} // dbt

#endif // #ifndef UTIL_H____
