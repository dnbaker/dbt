#ifndef CHARACTERHASH
#define CHARACTERHASH

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned int uint;

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <random>

using namespace std;



#if 0
class mersenneRNG {
public:
    mersenneRNG(uint32 maxval) : mtr(),n(maxval) {};
    uint32 operator()() {
        return mtr.randInt(n);
    }
    void seed(uint32 seedval) {
        mtr.seed(seedval);
    }
    void seed() {
        mtr.seed();
    }
    uint32 rand_max() {
        return n;
    }
private:
    MTRand mtr;
    int n;
};
#else
template<typename IType> struct MersenneWrapper;
template<> struct MersenneWrapper<uint32_t> {
    using type = std::mt19937;
};
template<> struct MersenneWrapper<uint64_t> {
    using type = std::mt19937;
};
#endif

#ifndef CONSTIF
#if __cplusplus >= 201703L
#define CONSTIF if constexpr
#else
#define CONSTIF if
#endif
#endif // #ifndef CONSTIF

template <typename hashvaluetype>
#if __cplusplus >= 201402L
constexpr
#endif
hashvaluetype maskfnc(int bits) {
    assert(bits>0);
    assert(bits<=sizeof(hashvaluetype)*8);
    hashvaluetype x = static_cast<hashvaluetype>(1) << (bits - 1);
    return x ^ (x - 1);
}

template <typename hashvaluetype = uint32, typename chartype = unsigned char, typename RNG=typename MersenneWrapper<std::make_unsigned_t<hashvaluetype>>::type>
class CharacterHash {
public:
    using mersenneRNG = RNG;
    CharacterHash(hashvaluetype maxval) {
        mersenneRNG randomgenerator(sizeof(hashvaluetype) == 4 ? maxval: (maxval>>32));
        mersenneRNG randomgeneratorbase(maxval>>32 ? maxval:0xFFFFFFFFU);
        for(size_t k =0; k<nbrofchars; ++k) {
            hashvalues[k] = static_cast<hashvaluetype>(randomgenerator());
            CONSTIF(sizeof(hashvaluetype) == 8)
                hashvalues[k] |= static_cast<hashvaluetype>(randomgeneratorbase()) << 32;
        }
    }

    CharacterHash(hashvaluetype maxval, uint32 seed1, uint32 seed2) {
        if(sizeof(hashvaluetype) <=4) {
            mersenneRNG randomgenerator(maxval);
            randomgenerator.seed(seed1);
            for(size_t k =0; k<nbrofchars; ++k)
                hashvalues[k] = static_cast<hashvaluetype>(randomgenerator());
        } else if (sizeof(hashvaluetype) == 8) {
            mersenneRNG randomgenerator(maxval>>32);
            mersenneRNG randomgeneratorbase((maxval>>32) ==0 ? maxval : 0xFFFFFFFFU);
            randomgenerator.seed(seed1);
            randomgeneratorbase.seed(seed2);
            for(size_t k =0; k<nbrofchars; ++k)
                hashvalues[k] = static_cast<hashvaluetype>(randomgeneratorbase())
                                | (static_cast<hashvaluetype>(randomgenerator()) << 32);
        } else throw runtime_error("unsupported hash value type");
    }

    enum {nbrofchars = 1 << ( sizeof(chartype)*8 )};

    hashvaluetype hashvalues[1 << ( sizeof(chartype)*8 )];
};

#endif
