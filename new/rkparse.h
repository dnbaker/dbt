#ifndef RK_PARSE_H__
#define RK_PARSE_H__
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <cstdlib>
#include "klib/khash.h"

namespace bwt {


namespace detail {


#ifndef PRIME_NUMBER
#define PRIME_NUMBER 27162335252586509
#endif

static constexpr uint64_t prime = UINT64_C(PRIME_NUMBER);


} // namespace detail
// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
static inline uint64_t kr_hash(const char *s) {
// TODO: consider unrolling.
    assert(*s);
    uint64_t hash = *s++;
    while(*s) hash = ((hash << 8) | *it++) % detail::PRIME;
    return hash;
}

static inline uint64_t kr_hash(const std::string &s) {return kr_hash(s.data());}
using std::size_t;
// class to maintain a window in a string and its KR fingerprint
template<bool pow2_wsize=false, uint64_t prime=1999999973ull>
struct KR_window_base {
  uint8_t *window;
  uint64_t hash;
  uint64_t tot_char;
  const uint64_t asize_pot;   // asize^(wsize-1) mod prime 
  int wsize;
  
  static constexpr int      asize = 256; // So that the multiplies can be expressed as bit shifts
  static constexpr uint64_t PRIME = prime;
  KR_window_base(int w): window(new uint8_t[w]), asize_pot(modexp(asize, w - 1, prime)), wsize(w) {
    if(window == nullptr) throw std::bad_alloc();
    if(pow2_wsize && w && (w & (w-1))) throw std::runtime_error("Window size must be a power of two if pow2_wsize is true!");
    
    // alloc and clear window
    reset();     
  }
  int       *begin()       {return window;}
  const int *begin() const {return window;}
  int       *end()       {return window + wsize;}
  const int *end() const {return window + wsize;}
  
  // init window, hash, and tot_char 
  void reset() {
    std::memset(window, 0, wsize * sizeof(*window));
    // init hash value and related values
    hash=tot_char=0;    
  }
  
  uint64_t addchar(char c) {
    if(pow2_wsize) {
        int k = tot_char++ & (wsize - 1);
    } else {
        int k = tot_char++ % wsize;
    }
        // complex expression to avoid negative numbers 
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution  
    hash = (asize*hash + c) % prime;      //  add char i 
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash; 
  }
  // debug only 
  string get_window() const {
    string w;
    for(int k = (tot_char-1) % wsize, i=k+1;i<k+1+wsize;w.append(1,window[i++%wsize]));
    return w;
  }
  
  ~KR_window_base() {delete[] window;}

};
using KR_window = KR_window_base<>;

struct word_stats {
  char *str;
  occ_int_t occ;
  word_int_t rank=0;
};
KHASH_MAP_INIT_INT64(stats, word_stats)

template<typename IType=std::uint64_t, bool pow2_wsize=false, typename=typename std::enable_if<std::is_integral<IType>::value>::type>
class lz77_t {
    khash_t(stats) h_;
    KR_window_base<pow2_wsize>  kr_;
    lz77_t(int wsize, size_t reserve=0): kr_(wsize) {
        std::memset(&h_, 0, sizeof(h_));
        if(reserve) kh_resize(stats, &h_, reserve);
    }
    ~lz77_t() {std::free(h_.keys); std::free(h_.vals); std::free(h_.flags);}
};

} // namespace bwt

#endif /* RK_PARSE_H__ */
