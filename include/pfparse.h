#ifndef PREFIX_FREE_PARSE_H__
#define PREFIX_FREE_PARSE_H__
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include "klib/khash.h"
#include "rollinghashcpp/rabinkarphash.h"

namespace dbt {

template <typename hashvaluetype = uint32, typename chartype =  unsigned char, unsigned wordsize=64>
class KarpRabinHashBits {
    // Modified from https://github.com/lemire/rollinghashcpp
    // The key difference is that wordsize is now templated
    // And the masking is only performed if nbits != the number of bits in the type
public:
    // myn is the length of the sequences, e.g., 3 means that you want to hash sequences of 3 characters
    // mywordsize is the number of bits you which to receive as hash values, e.g., 19 means that the hash values are 19-bit integers
    KarpRabinHashBits(int myn):  hashvalue(0), n(myn),
        hasher( maskfnc<hashvaluetype>(wordsize)),
        HASHMASK(maskfnc<hashvaluetype>(wordsize)),BtoN(1) {
        for (int i=0; i < n ; ++i) {
            BtoN *= B;
            if constexpr(wordsize != (CHAR_BIT * sizeof(hashvaluetype))) BtoN &= HASHMASK;
        }
    }

    // prepare to process a new string, you will need to call "eat" again
    void reset() {
      hashvalue = 0;
    }
    static constexpr bool is_full_word() {
        return wordsize == (CHAR_BIT * sizeof(hashvaluetype));
    }

    // this is a convenience function, use eat,update and .hashvalue to use as a rolling hash function
    template<class container>
    hashvaluetype  hash(container & c) {
        hashvaluetype answer(0);
        for(uint k = 0; k<c.size(); ++k) {
            hashvaluetype x(1);
            for(uint j = 0; j< c.size()-1-k; ++j) {
                x= (x * B);
                if constexpr(!is_full_word()) x &= HASHMASK;
            }
            x= (x * hasher.hashvalues[c[k]]);
            if constexpr(!is_full_word()) x &= HASHMASK;
            answer=(answer+x);
            if constexpr(!is_full_word()) answer &= HASHMASK;
        }
        return answer;
    }

    // add inchar as an input, this is used typically only at the start
    // the hash value is updated to that of a longer string (one where inchar was appended)
    void eat(chartype inchar) {
        hashvalue = (B*hashvalue +  hasher.hashvalues[inchar] );
        if constexpr(!is_full_word()) hashvalue &= HASHMASK;
    }

    // add inchar as an input and remove outchar, the hashvalue is updated
    // this function can be used to update the hash value from the hash value of [outchar]ABC to the hash value of ABC[inchar]
    void update(chartype outchar, chartype inchar) {
        hashvalue = (B*hashvalue +  hasher.hashvalues[inchar] - BtoN *  hasher.hashvalues[outchar]);
        if(!is_full_word()) hashvalue &= HASHMASK;
    }


    hashvaluetype hashvalue;
    int n;
    CharacterHash<hashvaluetype,chartype> hasher;
    const hashvaluetype HASHMASK;
    hashvaluetype BtoN;
    static constexpr hashvaluetype B=37;
};

KHASH_MAP_INIT_INT64(m, const char *)

class khmap {
    khash_t(m) map_;
public:
    khmap() {
        std::memset(&map_, 0, sizeof(map_));
    }
    void insert(uint64_t v, const char *s, size_t nelem) {
        khiter_t ki = kh_get(m, &map_, v);
        if(__builtin_expect(ki != kh_end(&map_), 0)) {
            throw std::runtime_error("Hash collision. Abort!");
        }
        int khr;
        ki = kh_put(m, &map_, v, &khr);
        char *s2 = static_cast<char *>(std::malloc(nelem + 1));
        if(!s2) throw std::bad_alloc();
        std::memcpy(s2, s, nelem);
        s2[nelem] = '\0';
        //auto ptr = &map_;
        kh_val(&map_, ki) = const_cast<const char *>(s2);
    }
    void insert(uint64_t v, const char *s) {insert(v, s, std::strlen(s));}
    khmap(khmap &&m) {std::memset(this, 0, sizeof(*this)); *this = std::move(m);}
    khmap &operator=(khmap &&m) {
        std::free(map_.keys);
        std::free(map_.vals);
        std::free(map_.flags);
        std::memcpy(this, &m, sizeof(m));
        m.map_ = khash_t(m){0,0,0,0,0,0,0};
        //std::memset(&m, 0, sizeof(m));
        return *this;
    }
    khmap(const khmap &) = delete;
    khmap& operator=(const khmap &) = delete;
    ~khmap() {
        for(khiter_t ki = 0; ki < map_.n_buckets; ++ki)
            if(kh_exist(&map_, ki))
                std::free(const_cast<char *>(map_.vals[ki]));
        std::free(map_.keys);
        std::free(map_.vals);
        std::free(map_.flags);
    }
};

using Hasher = KarpRabinHashBits<uint64_t>;

class HashPass {
    Hasher h_;
    khmap map_;
    uint64_t i_;
    std::deque<unsigned char> cstr_;
    // Prime 1?
public:
    static constexpr uint64_t PRIME = (1ull << 63) - 1;
    HashPass(unsigned wsz): h_(wsz), i_(0) {
    }
};


} // namespace dbt

#endif /* #ifndef PREFIX_FREE_PARSE_H__ */
