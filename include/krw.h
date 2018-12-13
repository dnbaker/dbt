#ifndef KR_WINDOW_H__
#define KR_WINDOW_H__
#include <cstdint>
#include <cstdlib>
#include <string>

namespace krw {
using std::uint8_t;
using std::uint64_t;
static constexpr uint64_t KRW_PRIME = 1999999973;
static constexpr int      KRW_ASIZE = 256;
static constexpr uint64_t make_asize_pot(int wsz) {
    uint64_t asize_pot = 1;
    for(int i=1;i<wsz;asize_pot <<= 8, asize_pot %= KRW_PRIME, ++i); // ugly linear-time power algorithm, but only done at construction.
    return asize_pot;
}


// class to maintain a window in a string and its KR fingerprint
struct KRWindow {
  // Based on the KR_window from Manzini's Big-BWT project: https://gitlab.com/manzai/Big-BWT / https://gitlab.com/manzai/Big-BWT/blob/master/newscan.cpp
  int wsize;
  int k_;
  std::string window;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime

  KRWindow(int w): wsize(w), k_(0), window(w, 0),
                   hash(0), tot_char(0), asize_pot(make_asize_pot(w))
  {
    // alloc and clear window
    reset();
  }
  KRWindow(const KRWindow&) = delete;
  KRWindow(KRWindow &&o): wsize(o.wsize), k_(o.k_), window(std::move(o.window)),
                          hash(o.hash), tot_char(o.tot_char), asize_pot(o.asize_pot)
  {
  }

  uint8_t       *begin() {
        return reinterpret_cast<uint8_t *>(window.data());
  }
  const uint8_t  *begin() const {
        return reinterpret_cast<const uint8_t *>(window.data());
  }

  uint8_t       *end()         {return reinterpret_cast<uint8_t *>(window.data() + window.size());}
  const uint8_t *end()   const {return reinterpret_cast<const uint8_t *>(window.data() + window.size());}

  // init window, hash, and tot_char
  void reset() {
    std::memset(window.data(), 0, window.size());
    // init hash value and related values
    hash=tot_char=k_=0;
  }

  uint64_t addchar(int c) {
    ++tot_char;
    if(++k_ == wsize)
        k_ = 0;
    // complex expression to avoid negative numbers
    hash += (KRW_PRIME - (window[k_]*asize_pot) % KRW_PRIME); // remove window[k] contribution
    window[k_] = c;
    hash = ((hash<<8) | c) % KRW_PRIME;      //  add char i
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash;
  }
  // debug only
  std::string get_window() const {
    std::string w;
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }

};

} // namespace krw

#endif
