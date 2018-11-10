#pragma once
#include "klib/khash.h"
#include "rollinghashcpp/rabinkarphash.h"

namespace dbt {

KHASH_MAP_INIT_INT64(64, uint64_t)
KHASH_MAP_INIT_INT(32, uint32_t)

class BlockTree {
    // Two passes per level:
    // 1. 
    uint64_t l_;
    unsigned r_; // length of each block
    unsigned w_; // window size
    khash_t(64) h_;
    std::FILE *ip_;
    const char *s_;
    size_t strind_;

    struct BlkData {
        uint64_t i_; // index
        uint64_t hash_;
    };
    static constexpr unsigned char SENTINEL = 0;
    KarpRabinHash<uint64_t> hasher_;
public:
#if 0
template <typename hashvaluetype = uint32, typename chartype =  unsigned char>
class KarpRabinHash {

public:
    // myn is the length of the sequences, e.g., 3 means that you want to hash sequences of 3 characters
    // mywordsize is the number of bits you which to receive as hash values, e.g., 19 means that the hash values are 19-bit integers
    KarpRabinHash(int myn, int mywordsize=19) :  hashvalue(0),n(myn),
        wordsize(mywordsize),
        hasher( maskfnc<hashvaluetype>(wordsize)),
        HASHMASK(maskfnc<hashvaluetype>(wordsize)),BtoN(1) {
#endif
    BlockTree(int wsz, uint64_t len, unsigned r): l_(len), r_(r), w_(wsz), ip_(nullptr), s_(0), strind_(0), hasher_(wsz, 64) {
        std::memset(&h_, 0, sizeof(h_));
    }
    std::vector<BlkData> first_pass(unsigned lvl) {
        std::vector<BlkData> ret;
        unsigned blen = len / r_;
        for(unsigned k = lvl; --k; blen /= r_);
        int c;
        unsigned nfilled;
        if(l_ > w_) throw std::runtime_error("Can't hash a sequence shorter than the window.");
        for(nfilled = 0; nfilled < w_ && c = getnext() != EOF; hasher.eat(c));
        ret.emplace_back(hash_);
        nfilled = 0;
        while((c = this->getnext()) != EOF) {
            hasher_.eat(c);
            if(++nfilled == w_) {
                ret.push_back(hasher_.hashvalue);
                nfilled = 0;
            }
        }
        while(nfilled++ < w_)
            hasher_.eat(SENTINEL);
        ret.push_back(hasher_.hashvalue);
        return ret;
    }
    void rewind() {
        if(ip_) std::rewind(ip_);
        else         strind_ = 0;
    }
    unsigned char getnext() noexcept {
        if(ip_) {
            return std::fgetc(ip_);
        } else {
            if(__builtin_expect(strind_ < l_, 1))
                return s_[strind_++];
            else return EOF;
        }
    }
};

} // namespace dbt
