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

namespace dbt {


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
    void free_map() {
        std::free(map_.keys);
        std::free(map_.flags);
        std::free(map_.vals);
        std::memset(&map_, 0, sizeof(map_));
    }
    void free_values() {
        for(khiter_t ki = 0; ki < map_.n_buckets; ++ki)
            if(kh_exist(&map_, ki))
                std::free(const_cast<char *>(map_.vals[ki].s_));
    }
    void free() {
        free_values();
        free_map();
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
    const khash_t(m) *map() const {return &map_;}
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
    ~khmap() {this->free();}
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
                    kh_val(&map_, ki).s_ = util::strdup(kh_val(&map_, ki).s_);
                }
            }
        }
        return *this;
    }
};

enum {Dollar = 0xFF, EndOfWord = 0x00};

using Hasher = KarpRabinHashBits<uint64_t>;

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
    struct StrCmp {
        bool constexpr operator()(const char *s, const char *s2) const {
            return std::strcmp(s, s2) < 0;
        }
    };

    std::set<const char *, StrCmp> pointers;
    std::string dictpath = prefix; dictpath += DEXT;
    std::string occpath = prefix; occpath   += OCCEXT;
    util::FpWrapper<PointerType> dfp(dictpath.data(), "wb");
    util::FpWrapper<PointerType> ocfp(occpath.data(), "wb");
    for(auto &h: hp) {
        if(h.words) {
            tot += h.parse - (i++ ? h.w(): 0);
        }
        for(size_t j = 0; j < h.map_->capacity(); ++j) {
            if(kh_exist(h.map_->map(), j)) {
                auto ip = pointers.insert(kh_val(h.map_->map(), j).s_);
                if(ip.second == false) {
                    std::free(const_cast<char *>(kh_val(h.map_->map(), j).s_));
                    kh_val(h.map_->map(), j).s_ = nullptr;
                }
            }
        }
        *map |= *h.map_;
        h.map_->free_map();
        // Can these resources be freed now?
    }
    std::vector<const char *> dictset(pointers.begin(), pointers.end());
    pointers.clear();
    for(auto &h: hp) h.map_->free();
    for(auto &x: dictset)
        x = util::strdup(x);
    std::sort(dictset.begin(), dictset.end(), [](const char *a, const char *b) {return std::strcmp(a, b) < 0;});
    uint64_t totDWord = map->size();
    size_t wrank = 0;
    // Write dictionary occurrences
    for(const auto c: dictset) {
        size_t sl = std::strlen(c), s = dfp.write(c);
        auto hv = hp[0].hash(c);
        if(s != sl) throw std::runtime_error("Could not write to dict file.");
        dfp.write(EndOfWord);
        std::free(const_cast<char *>(c));
        hit_t &h = map->operator[](hv);
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
    util::FpWrapper<std::FILE *> pfp;
    pfp.open((std::string(prefix) + PEXT).data(), "wb");
    util::FpWrapper<std::FILE *> tfp;
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



namespace lut {

//print(", ".join(map(str, (1 if x in 'Aa' else 2 if x in 'Cc' else 4 if x in 'Gg' else 8 if x in 'Tt' else 0 if x in '\t\n ' else 15 for x in map(chr, range(256))))))
static constexpr uint8_t lut4b [] {
    15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 0, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1, 15, 2, 15, 15, 15, 4, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1, 15, 2, 15, 15, 15, 4, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

static constexpr uint8_t lut2b [] {
    //print(", ".join(map(str, (0 if x in 'Aa' else 1 if x in 'Cc' else 2 if x in 'Gg' else 3 if x in 'Tt' else 4 for x in map(chr, range(256))))))
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

static constexpr uint8_t printable [] {
    33, 33, 33, 33, 33, 33, 33, 33, 33, 9, 10, 11, 12, 13, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33
};

} // namespace lut

template<typename PointerType>
class HashPass {
    using FType = util::FpWrapper<PointerType>;
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
    auto hash(const C &c) const {return h_.hash(c);}

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
    void open_ifp(const char *path, const char *mode="rb") {
        ifp_.open(path, mode);
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
    int filternextchar(uint8_t *tbl=lut::printable) {
        int nc;
        while((nc = nextchar()) != EOF && (nc = tbl[nc]) != 0);
        return nc;
    }
    void fill(size_t nelem=0, size_t start=0) {
        skip = 0, parse = 0, words = 0;
        assert(ifp_.is_open());
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
        } else cstr_.push_front(Dollar);
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
