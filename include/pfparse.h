#ifndef PREFIX_FREE_PARSE_H__
#define PREFIX_FREE_PARSE_H__
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cstdio>
#include <deque>
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif
#include "klib/khash.h"
#include "rollinghashcpp/rabinkarphash.h"
#include <sys/stat.h>

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
            if(!is_full_word()) BtoN &= HASHMASK;
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

struct hit_t {
    const char *s_;
    uint32_t rank_;
    uint32_t  occ_;
};

KHASH_MAP_INIT_INT64(m, hit_t)

class khmap {
    khash_t(m) map_;
public:
    khmap() {
        std::memset(&map_, 0, sizeof(map_));
    }
    // TODO: consider compressing strings by storing in 4-bits per character
    void insert(uint64_t v, const char *s, size_t nelem) {
        khiter_t ki = kh_get(m, &map_, v);
        if(ki != kh_end(&map_)) {
            if(__builtin_expect(std::strcmp(s, map_.vals[ki].s_), 0))
                throw std::runtime_error("Hash collision. Abort!");
            if(__builtin_expect(++map_.vals[ki].occ_ == 0, 0))
                throw std::runtime_error("Overflow in occurrence count");
        } else {
            int khr;
            ki = kh_put(m, &map_, v, &khr);
            char *s2 = static_cast<char *>(std::malloc(nelem + 1));
            if(!s2) throw std::bad_alloc();
            std::memcpy(s2, s, nelem);
            s2[nelem] = '\0';
            kh_val(&map_, ki) = hit_t{const_cast<const char *>(s2), 0, 1};
        }
    }
    void insert(uint64_t v, const char *s) {insert(v, s, std::strlen(s));}
    template<typename T>
    void insert(uint64_t v, const T &s) {
        std::string t(s.begin(), s.end());
        insert(v, t.data(), t.size());
    }
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
                std::free(const_cast<char *>(map_.vals[ki].s_));
        std::free(map_.keys);
        std::free(map_.vals);
        std::free(map_.flags);
    }
};

enum {Dollar = 0xFF};

using Hasher = KarpRabinHashBits<uint64_t>;

template<typename PointerType>
class FpWrapper {
    PointerType ptr_;
    std::vector<char> buf_;
public:
    using type = PointerType;
    FpWrapper(type ptr=nullptr): ptr_(ptr), buf_(BUFSIZ) {}
    static constexpr bool is_gz() {
        return std::is_same_v<PointerType, gzFile>;
    }
    static constexpr bool maybe_seekable() {
        return !std::is_same_v<PointerType, gzFile>;
    }
    bool seekable() const {
        if constexpr(is_gz()) return false;
        struct stat s;
        ::fstat(::fileno(ptr_), &s);
        return !S_ISFIFO(s.st_mode);
    }
    auto write(const void *buf, size_t n) {
        if constexpr(is_gz()) {
            return gzwrite(ptr_, buf, n);
        } else return ::write(::fileno(ptr_), buf, n);
    }
    template<typename T>
    auto write(T val) {
        this->write(&val, sizeof(val));
    }
    auto resize_buffer(size_t newsz) {
        if constexpr(!is_gz()) {
            buf_.resize(newsz);
            std::setvbuf(ptr_, buf_.data(), buf_.size());
        } else {
            gzbuffer(ptr_, newsz);
        }
    }
    void close() {
        if constexpr(is_gz())
            gzclose(ptr_);
        else
            fclose(ptr_);
    }
    void open(const char *path, const char *mode) {
        if constexpr(is_gz()) {
            ptr_ = gzopen(path, mode);
        } else {
            ptr_ = fopen(path, mode);
        }
        if(ptr_ == nullptr)
            throw std::runtime_error(std::string("Could not open file at ") + path + " with mode" + mode);
    }
    bool is_open() const {return ptr_ != nullptr;}
    auto eof() {
        if constexpr(is_gz())
            return gzeof(ptr_);
        else
            std::feof(ptr_);
    }
    ~FpWrapper() {
        close();
    }
};

template<typename PointerType=std::FILE*>
class HashPass {
    using FType = FpWrapper<PointerType>;
    Hasher h_;
    uint64_t i_;
    std::deque<unsigned char> cstr_;
    std::deque<unsigned char> q_;
    FType safp_, lastfp_, pafp_, ifp_, hafp_;
    const char *prefix_;
    std::vector<char> buf_;
    size_t             bi_;
    // Prime 1?
public:
    khmap *map_;



    static constexpr uint64_t LARGE_PRIME = (1ull << 63) - 1;
    static constexpr uint64_t SMALL_PRIME = 109829;
    static constexpr uint64_t MEDIUMPRIME = 1999999973;
    static constexpr uint64_t WINDOW_MOD  = 100;


    HashPass(unsigned wsz, const char *pref="default_prefix", int compression=6): h_(wsz), i_(0), prefix_(pref), bi_(0) {
        std::string safn = pref; safn += ".sa";
        std::string mode = "wb";
        if(safp_.is_gz()) {
            if(compression == 0) mode = "wT";
            else mode += std::to_string(compression);
        }
        safp_.open(safn.data(), mode.data());
        safn = pref; safn += ".last";
        lastfp_.open(safn.data(), mode.data());
        safn = pref; safn += ".prs";
        pafp_.open(safn.data(), mode.data());
        safn = pref; safn += ".hv";
        hafp_.open(safn.data(), mode.data());
        buf_.resize(1 << 14);
    }
    int nextchar() {
        if(__builtin_expect(bi_ == buf_.size(), 0)) {
            size_t n = ifp_->read(buf_.data(), buf_.size());
            if(!n) return EOF;
            if(n != buf_.size())
                buf_.resize(n);
            bi_ = 0;
        }
        return buf_[bi_++];
    }
    void fill(size_t nelem=0, size_t start=0) {
        size_t skip = 0, parse = 0, words = 0;
        if(start) {
            int c;
            for(;(c = nextchar()) != EOF && q_.size() != h_.n;q_.push_back(c)) {
                h_.eat(c);
                q_.push_back(c);
                if(++skip == nelem + start) {std::fprintf(stderr, "Warning: sequence too short\n"); return;}
            }
            if(q_.size() != w()) {std::fprintf(stderr, "Warning: could not fill. Returning early\n"); return;}
            while(skip < w() || h_.hashvalue % WINDOW_MOD) {
                if((c = nextchar()) == EOF) {std::fprintf(stderr, "Warning: sequence too short2\n"); return;}
                if(++skip == nelem + start) {std::fprintf(stderr, "Warning: sequence too short3\n"); return;}
                h_.update(q_.front(), c);
                q_.pop_front();
                q_.push_back(c);
            }
            parse = h_.n;
            skip -= h_.n;
            assert(h_.n == q_.size());
            cstr_ = q_;
        }
        else
            cstr_.push_front(Dollar);
        size_t pos = start;
        if(pos) pos += skip + w();
        for(int c; (c = nextchar()) != EOF;) {
            ++parse;
            h_.update(q_.front(), c);
            q_.pop_front();
            q_.push_back(c);
            cstr_.push_back(c);
            if(h_.hashvalue % WINDOW_MOD == 0) {
                ++words;
                update(pos);
                if(skip + parse == nelem + w()) break;
            }
        }
        cstr_.insert(cstr_.end(), w(), Dollar);
        update(pos);
    }
    void update(uint64_t &pos) {
        auto hv = h_.hash(cstr_);
        hafp_.write(hv); 
        map_->insert(hv, cstr_);
        lastfp_.write(cstr_[cstr_.size() - 1 - w()]);
        if(pos==0) pos = cstr_.size()-1;
        else       pos += cstr_.size() - w();
        if(safp_.is_open())
            safp_.write(pos);
        cstr_.erase(cstr_.begin(), cstr_.size() - w());
    }
    auto w() const {return h_.n;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modlp(T v) {return v & LARGE_PRIME;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modsp(T v) {return v & MEDIUMPRIME;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modmp(T v) {return v & SMALL_PRIME;}
};


} // namespace dbt

#endif /* #ifndef PREFIX_FREE_PARSE_H__ */
