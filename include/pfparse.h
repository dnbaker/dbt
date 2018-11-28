#ifndef PREFIX_FREE_PARSE_H__
#define PREFIX_FREE_PARSE_H__
#include <string>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <deque>
#include <set>
#include "util.h"
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
        if constexpr(!is_full_word()) hashvalue &= HASHMASK;
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
    uint64_t rank_;
    uint32_t  occ_;
};

KHASH_MAP_INIT_INT64(m, hit_t)

class khmap {
    khash_t(m) map_;
public:
    khmap() {
        std::memset(&map_, 0, sizeof(map_));
    }
    void free() {
        for(khiter_t ki = 0; ki < map_.n_buckets; ++ki)
            if(kh_exist(&map_, ki))
                std::free(const_cast<char *>(map_.vals[ki].s_));
        std::free(map_.keys);
        std::free(map_.flags);
        std::free(map_.vals);
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
        this->free();
        std::memcpy(this, &m, sizeof(m));
        m.map_ = khash_t(m){0,0,0,0,0,0,0};
        //std::memset(&m, 0, sizeof(m));
        return *this;
    }
    hit_t &operator[](uint64_t v) {
        khiter_t ki = kh_get(m, &map_, v);
        if(ki == kh_end(&map_)) throw std::out_of_range(std::string("Missing key ") + std::to_string(v));
        return kh_val(&map_, ki);
    }
    const hit_t &operator[](uint64_t v) const {
        khiter_t ki = kh_get(m, &map_, v);
        if(ki == kh_end(&map_)) throw std::out_of_range(std::string("Missing key ") + std::to_string(v));
        return kh_val(&map_, ki);
    }
    khmap(const khmap &) = delete;
    khmap& operator=(const khmap &) = delete;
    ~khmap() {
        this->free();
    }
    size_t size() const {return kh_size(&map_);}
    size_t capacity() const {return map_.n_buckets;}
    khmap &operator|=(const khmap &o) {
        for(size_t i = 0; i < o.map_.n_buckets; ++i) {
            if(kh_exist(&o.map_, i)) {
                int khr;
                khiter_t ki = kh_put(m, &map_, kh_key(&o.map_, i), &khr);
                if(khr == 0) {
                    if(std::strcmp(kh_val(&o.map_, i).s_, kh_val(&map_, ki).s_)) throw std::runtime_error("ZOMG");
                    kh_val(&map_, ki).occ_ += kh_val(&o.map_, i).occ_;
                } else {
                    kh_val(&map_, ki) = kh_val(&o.map_, i);
                }
            }
        }
    }
};

enum {Dollar = 0xFF, EndOfWord = 0x00};

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
    template<typename T>
    auto read(T *val) {
        return this->read(val, sizeof(T));
    }
    auto read(void *ptr, size_t nb) {
        if constexpr(is_gz())
            return gzread(ptr_, ptr, nb);
        else
            return std::fread(ptr, 1, nb, ptr_);
    }
    auto write(const char *s) {
        if constexpr(is_gz())
            return gzputs(ptr_, s);
        else
            return std::fputs(s, ptr_);
    }
    auto resize_buffer(size_t newsz) {
        if constexpr(!is_gz()) {
            buf_.resize(newsz);
            std::setvbuf(ptr_, buf_.data(), buf_.size());
        } else {
            gzbuffer(ptr_, newsz);
        }
    }
    void seek(size_t pos) {
        if constexpr(is_gz())
            gzseek(ptr_, pos);
        else
            std::fseek(ptr_, pos, SEEK_SET);
    }
    void close() {
        if constexpr(is_gz())
            gzclose(ptr_);
        else
            fclose(ptr_);
        ptr_ = nullptr;
    }
    auto write(void *buf, size_t nelem) {
        if constexpr(is_gz())
            return gzwrite(ptr_, buf, nelem);
        else
            return std::fwrite(buf, 1, nelem, ptr_);
    }
    template<typename T>
    void write(T val) {this->write(&val, sizeof(val));}
    void open(const char *path, const char *mode) {
        if(ptr_) close();
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
        if(ptr_) close();
    }
};

struct StrCmp {
    bool constexpr operator()(const char *s, const char *s2) const {
        return std::strcmp(s, s2) < 0;
    }
};

char *strdup(const char *s) {
    size_t l = std::strlen(s) + 1;
    char *ret = static_cast<char *>(std::malloc(l));
    if(!ret) throw std::bad_alloc();
    std::memcpy(ret, s, l);
    return ret;
}

const std::string DEXT    = ".dict";
const std::string PEXT    = ".parse";
const std::string OCCEXT  = ".occ";

template<typename PointerType=std::FILE*> class HashPass;

template<typename PointerType=std::FILE*>
void merge_hashpasses(const char *prefix, const std::vector<HashPass<PointerType>> &hp, khmap *map) {
/*
    When this completes, we have a dictionary, occurrence, and a parse file.
    The suffix array pieces have been generated for each subset, but not yet concatenated.
    This is for merging of multiple efforts on a single genome.

    dictionary: .dict
    occ:        .occ
    parse:      .parse
    sa (unmerged) .<int>.sa // Can be concatenated
    tot:       (in memory)  // Total number of characters in input, minus the w extra

    Merging multiples of these will essentially require:
        1. Getting lengths of full corpora (value of 'tot' here)
        2. Assigning indexes to each genome and offsets to increment by
        3. Incrementing all position information by each respective offset to get to the global index
        4. Doing this for the sa files before concatenating
 */
    assert(map);
    size_t tot = 0, i = 0;
    std::set<const char *, StrCmp> pointers;
    std::string dictpath = prefix; dictpath += DEXT;
    std::string occpath = prefix; occpath   += OCCEXT;
    FpWrapper dfp(dictpath.data(), "wb");
    FpWrapper ocfp(occpath.data(), "wb");
    for(const auto &h: hp) {
        if(h.words) {
            tot += h.parse - (i++ ? h.w(): 0);
        }
        for(size_t j = 0; j < h.map_.n_buckets; ++j) {
            if(kh_exist(&h.map_, j))
                pointers.insert(kh_val(&h.map_, j).s_);
        }
        *map |= *h.map_;
        h.map_->free();
        // Can these resources be freed now?
    }
    std::vector<char *> dictset(pointers.begin(), pointers.end());
    pointers.clear();
    for(auto &x: dictset) {
        x = strdup(x);
    }
    std::sort(dictset.begin(), dictset.end(), [](const char *a, const char *b) {return std::strcmp(a, b) < 0;});
    uint64_t totDWord = kh_size(map);
    size_t wrank = 0;
    // Write dictionary occurrences
    for(const auto c: dictset) {
        size_t sl = std::strlen(c), s = dfp.write(c);
        auto hv = hp[0].hash(c);
        if(s != sl) throw std::runtime_error("Could not write to dict file.");
        dfp.write(EndOfWord);
        std::free(c);
        hit_t &h = map[hv];
        assert(h.occ_);
        ocfp.write(h.occ_);
        assert(h.rank_ == 0);
        h.rank_ = ++wrank;
    }
    dfp.close();
    ocfp.close();
    size_t chunknum = 0;
#if MAKE_HISTOGRAM
    std::vector<uint64_t> ranks(map.size() + 1);
#endif
    FpWrapper<std::FILE *> pfp;
    pfp.open((std::string(prefix) + PEXT).data(), "wb");
    FpWrapper<std::FILE *> tfp;
    tfp.open((std::string(prefix) + '.' + std::to_string(chunknum) + PEXT).data(), "rb");
    uint64_t hv;
    FOREVER {
        if(tfp.read(&hv) != sizeof(hv)) {
            if(chunknum + 1 == hp.size())
                break;
            tfp.open((std::string(prefix) + '.' + std::to_string(++chunknum) + PEXT).data(), "rb");
            continue;
        }
        uint64_t rank = map->operator[](hv).rank_;
#if MAKE_HISTOGRAM
        ++ranks[rank];
#endif
        pfp.write(rank);
    }
}

template<typename PointerType>
class HashPass {
    using FType = FpWrapper<PointerType>;
public:
    Hasher h_;
private:
    uint64_t i_;
    std::deque<unsigned char> cstr_;
    std::deque<unsigned char> q_;
    FType safp_, lastfp_, pafp_, ifp_;
    const char *prefix_;
    std::vector<char> buf_;
    size_t             bi_;
    // Prime 1?
public:
    size_t skip, parse, words;
    khmap *map_;

    template<typename C>
    auto hash(const C &c) {return h_.hash(c);}

    static constexpr uint64_t LARGE_PRIME = (1ull << 63) - 1;
    static constexpr uint64_t SMALL_PRIME = 109829;
    static constexpr uint64_t MEDIUMPRIME = 1999999973;
    static constexpr uint64_t WINDOW_MOD  = 100;

    HashPass(unsigned wsz, size_t nchunks=1, size_t chunknum=0, const char *pref="default_prefix", int compression=6): h_(wsz), i_(0), prefix_(pref), bi_(0) {
        std::string safn = pref;
        if(nchunks > 1) safn += '.', safn += std::to_string(chunknum);
        safn += ".sa";
        std::string mode = "wb";
        if(safp_.is_gz()) {
            if(compression == 0) mode = "wT";
            else mode += std::to_string(compression);
        }
        safp_.open(safn.data(), mode.data());
        safn = pref;
        if(nchunks > 1) safn += '.', safn += std::to_string(chunknum);
        safn += ".last";
        lastfp_.open(safn.data(), mode.data());
        safn = pref;
        if(nchunks > 1) safn += '.', safn += std::to_string(chunknum);
        safn += ".parse";
        pafp_.open(safn.data(), mode.data());
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
        skip = 0, parse = 0, words = 0;
        if(start) {
            ifp_.seek(start);
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
        pafp_.write(hv);
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
    static constexpr auto modmp(T v) {return v & MEDIUMPRIME;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modsp(T v) {return v & SMALL_PRIME;}
};


} // namespace dbt

#endif /* #ifndef PREFIX_FREE_PARSE_H__ */
