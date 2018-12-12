#ifndef PREFIX_FREE_PARSE_H__
#define PREFIX_FREE_PARSE_H__
#include <algorithm>
#include <cstdio>
#include <deque>
#include <functional>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#include "util.h"
#include "klib/khash.h"
#include "logutil.h"
template<typename T> class TD;

#ifndef INLINE
#  if __GNUC__ || __clang__ || INLINE_OVERRIDE
#    define INLINE __attribute__((always_inline))
#  else
#    define INLINE
#  endif
#endif

#include "include/rabinkarphash.h"

namespace dbt {

namespace {
template<typename T>
static std::string stringify(const T &x) {
    auto it = std::begin(x);
    std::string ret = std::to_string(*it);
    while(it != std::end(x))
        ret += ',', ret += std::to_string(*it++);
    return ret;
}
} // anonymous namespace

struct hit_t {
    const char *s_;
    uint64_t rank_;
    uint32_t  occ_;
    operator const char *() const {return s_;}
};

KHASH_MAP_INIT_INT64(m, hit_t)

template<typename T>
void print_as_values(const T &cstr_) {
    std::string tmp;
    for(const auto c: cstr_) tmp += std::to_string(int(c)), tmp += ',';
    tmp.back() = '\n';
    LOG_DEBUG("cstr values: %s\n", tmp.data());
}

class khmap {
    khash_t(m) map_;
public:
    using value_type     = std::decay_t<decltype(*map_.vals)>;
    using reference_type = value_type &;
    using pointer_type   = value_type *;
    using key_type       = std::decay_t<decltype(*map_.keys)>;

    khmap() {std::memset(this, 0, sizeof(map_));}
    void free_map() {
        std::free(map_.keys);
        std::free(map_.flags);
        std::free(map_.vals);
        std::memset(this, 0, sizeof(map_));
    }
    void free() {
        free_values();
        free_map();
    }
    auto get(uint64_t v) const {
        return kh_get(m, _h(), v);
    }
    auto       &val(uint64_t ki)       {return kh_val(_h(), ki);}
    const auto &val(uint64_t ki) const {return kh_val(_h(), ki);}
    key_type key(uint64_t ki) const {return kh_key(_h(), ki);}
    hit_t &operator[](uint64_t v) {
#define HIT_CORE\
        if(auto ki = get(v); unlikely(is_end(v)))\
            throw std::out_of_range(std::string("Missing key ") + std::to_string(v));\
        else\
            return val(ki);
        HIT_CORE
    }
    const hit_t &operator[](uint64_t v) const {
        HIT_CORE
#undef HIT_CORE
    }
    khmap(const khmap &) = delete;
    khmap& operator=(const khmap &) = delete;
    ~khmap() {
#if !NDEBUG
        std::fprintf(stderr, "Destroying khmap at %p\n", static_cast<void *>(this));
#endif
        this->free();
    }
    size_t size() const {return kh_size(&map_);}
    size_t capacity() const {return map_.n_buckets;}
    bool is_end(khint_t ki) const {return ki == capacity();}
    bool exist(khint_t ki) const {return kh_exist(_h(), ki);}
    void print_all() const {
        for_each([](auto x, const hit_t &y) {
            std::fprintf(stderr, "String %s has hash %lu\n", y.s_, x);
        });
    }
    size_t memory_usage(bool add_strlens=false) const {
        size_t ret = sizeof(*this) + size() * (sizeof(value_type) + sizeof(key_type)) +
        __ac_fsize(map_.n_buckets) * sizeof(uint32_t);
        if(add_strlens)
            for_each_val([&](const hit_t &h) {std::fprintf(stderr, "h.s: %s/%zu\n", h.s_, std::strlen(h.s_)); ret += std::strlen(h.s_);});
        return ret;
    }
    template<class T, class GetterFunc, class BinaryOperation=std::plus<T>>
    T accumulate(khint_t start, khint_t end, T init, 
                 BinaryOperation op, const GetterFunc &gfunc) const
    {
        //end = end ? end: capacity();
        assert(start < end);
        while(start < end) {
            if(exist(start))
                init = op(std::move(init), gfunc(start++)); // std::move since C++20
        }
        return init;
    }
    template<class T, class BinaryOperation=std::plus<T>>
    INLINE T vacc(khint_t start, khint_t end, T init=T(), BinaryOperation op=BinaryOperation{}) const {
        return accumulate(start, end, init, op, [this](khint_t x) {return val(x);});
    }
    template<class T, class BinaryOperation=std::plus<T>>
    INLINE T kacc(khint_t start, khint_t end, T init=T(), BinaryOperation op=BinaryOperation{}) const {
        return accumulate(start, end, init, op, [this](khint_t x) {return val(x);});
    }
    template<class T, class BinaryOperation=std::plus<T>>
    INLINE T kvacc(khint_t start, khint_t end, T init, BinaryOperation op) const {
        return accumulate(start, end, init, op, [this](khint_t x) {return std::make_pair(key(x), std::ref(val(x)));});
    }
    auto put(uint64_t k) {
        int khr;
        return kh_put(m, _h(), k, &khr);
    }
    bool contains_key(uint64_t k) const {return !is_end(get(k));}
    // TODO: consider compressing strings by storing in 4-bits per character
    void insert(uint64_t v, const char *s, size_t nelem) {
#if !NDEBUG
        print_all();
#endif
        if(auto ki = get(v); !is_end(ki)) {
            assert(map_.vals[ki].s_);
            if(unlikely(std::strcmp(s, map_.vals[ki].s_)))
                throw std::runtime_error("Hash collision. Abort!");
            if(unlikely(++map_.vals[ki].occ_ == 0)) // This line both increments and checks! Do not remove.
                throw std::runtime_error("Overflow in occurrence count");
        } else {
#if !NDEBUG
            std::fprintf(stderr, "String: %s/%zu.\n", s, nelem);
#endif
            ki = put(v);
            char *s2 = static_cast<char *>(std::malloc(nelem + 1));
            if(!s2) throw std::bad_alloc();
            std::memcpy(s2, s, nelem);
            s2[nelem] = '\0';
            val(ki) = hit_t{const_cast<const char *>(s2), 0, 1};
            std::set<char> set(s2, s2 + nelem);
#if 0
            assert(std::accumulate(s2, s2 + nelem, true, [](bool v, const signed char c) -> bool {
                return v && ((c >= 0) || c == util::Dollar);
            }) ||
                  !std::fprintf(stderr, "Contents of set: %s\n", stringify(set).data()));
#endif
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
    operator       khash_t(m) *()        {return reinterpret_cast<khash_t(m) *>(this);}
    operator const khash_t(m) *()  const {return reinterpret_cast<const khash_t(m) *>(this);}
    khash_t(m) *_h() {return reinterpret_cast<khash_t(m) *>(this);}
    const khash_t(m) *_h() const {return reinterpret_cast<const khash_t(m) *>(this);}
        
    khmap &operator=(khmap &&m) {
        this->free();
        std::memcpy(this, &m, sizeof(m));
        m.map_ = khash_t(m){0,0,0,0,0,0,0};
        //std::memset(&m, 0, sizeof(m));
        return *this;
    }
    void swap(khmap &m) {
        //char buf[sizeof(m)];
        std::swap_ranges(reinterpret_cast<uint8_t *>(this), reinterpret_cast<uint8_t *>(this + sizeof(*this)), reinterpret_cast<uint8_t *>(std::addressof(m)));
        //std::memcpy(buf, &m, sizeof(m));
        //std::memcpy(&m, this, sizeof(m));
        //std::memcpy(this, buf, sizeof(*this));
    }
    void assert_nonnull() const {
#if 0
        for_each([](auto k, const hit_t &h) {assert(h.s_);});
#endif
    }
    template<typename Func>
    void __for_each_range(const Func &func) const {
        for(khiter_t ki = 0; ki < capacity(); ++ki) if(exist(ki)) func(ki);
    }
    template<typename Func>
    void __for_each_range(const Func &func) {
        for(khiter_t ki = 0; ki < capacity(); ++ki) if(exist(ki)) func(ki);
    }
    void free_values() {
        for_each_val([](hit_t &h) {std::free(const_cast<char *>(h.s_)); h.s_ = nullptr;});
    }
    template<typename Func>
    void __for_each(const Func &func) {
        for(khiter_t ki = 0; ki < capacity(); ++ki) if(exist(ki)) func(ki);
    }
    template<typename Func>
    void for_each_key(const Func &func) {__for_each_range([&](auto ki) {return func(this->key(ki));});}
    template<typename Func>
    void for_each_val(const Func &func) {__for_each_range([&](auto ki) {return func(this->val(ki));});}
    template<typename Func>
    void for_each_key(const Func &func) const {__for_each_range([&](auto ki) {return func(this->key(ki));});}
    template<typename Func>
    void for_each_val(const Func &func) const {__for_each_range([&](auto ki) {return func(this->val(ki));});}
    template<typename Func>
    void for_each(const Func &func) const {
        for(khiter_t ki = 0; ki < kh_end(&map_); ++ki)
            if(exist(ki))
                func(key(ki), val(ki));
    }
    template<typename Func>
    void for_each(const Func &func) {
        for(khiter_t ki = 0; ki < kh_end(&map_); ++ki)
            if(exist(ki))
                func(key(ki), val(ki));
    }
    khmap &operator|=(const khmap &o) {
        assert_nonnull();
        o.assert_nonnull();
        o.for_each([&](const auto key, const hit_t &h) {
            if(auto ki = get(key); is_end(ki)) {
                ki = put(key);
                val(ki) = h;
                val(ki).s_ = util::strdup(h.s_);
            } else val(ki).occ_ += h.occ_;
        });
        return *this;
    }
};


using Hasher = KarpRabinHashBits<uint64_t>;

const std::string DEXT    = ".dict";
const std::string PEXT    = ".parse";
const std::string OCCEXT  = ".occ";

template<typename PointerType=std::FILE*> class HashPass;

template<typename PointerType=std::FILE*>
void merge_hashpasses(const char *prefix, const std::vector<HashPass<PointerType>> &hp, const std::vector<std::string> &paths, khmap *map) {
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

    std::string dictpath = prefix; dictpath += DEXT;
    std::string occpath = prefix; occpath   += OCCEXT;
    util::FpWrapper<PointerType> dfp(dictpath.data(), "wb");
    util::FpWrapper<PointerType> ocfp(occpath.data(), "wb");
    if(hp.size() == 1) {
        tot = hp[0].parse;
        map->swap(*hp[0].map_);
        LOG_DEBUG("Swapping instead of merging\n");
    } else {
        LOG_DEBUG("Merging instead of swapping with %zu hps\n", hp.size());
        for(auto &h: hp) {
            if(h.words) {
                tot += h.parse - (i++ ? h.w(): 0);
            }
            *map |= *h.map_;
            //h.map_->free();
            // Can these resources be freed now?
        }
    }
    const char **dictset = static_cast<const char **>(
        std::malloc(map->size() * sizeof(const char *))), **it = dictset;
    map->for_each([&](auto k, const auto &h) {assert(h.s_);*it++ = h.s_;});
    std::sort(dictset, dictset + map->size(), [](const char *a, const char *b) {return std::strcmp(a, b) < 0;});
    uint64_t totDWord = map->size();
    LOG_DEBUG("totDWord: %" PRIu64 "\n", totDWord);
    uint64_t wrank = 0;
    // Write dictionary occurrences
    for(it = dictset; it < dictset + map->size();) {
        const char *c = *it++;
        size_t sl = std::strlen(c), s = dfp.write(c);
        auto hv = hp[0].hash(c);
        //std::fprintf(stderr, "string %s/%zu has has %" PRIu64 " hash value\n", c, sl, hv);
        if(s != sl) {
            char buf[256];
            std::sprintf(buf, "Could not write to dict file. string: %s s: %zu. sl: %zu", c, s, sl);
            throw std::runtime_error(buf);
        }
        dfp.write(util::EndOfWord);
        hit_t &h = map->operator[](hv);
        assert(h.occ_);
        ocfp.write(h.occ_);
        assert(h.rank_ == 0);
        h.rank_ = ++wrank;
    }
    std::free(dictset);
    dfp.close();
    ocfp.close();
    size_t ranknum = 0;
    util::FpWrapper<std::FILE *> pfp;
    pfp.open((std::string(prefix) + PEXT).data(), "wb");
    util::FpWrapper<std::FILE *> tfp;
    LOG_DEBUG("Paths! %zu of them\n", paths.size());
#if MAKE_HISTOGRAM
    std::vector<uint64_t> ranks(map.size() + 1);
#endif
    for(const auto &pref: paths) {
        // std::fprintf(stderr, "Opening file at %s\n", (pref + PEXT).data());
        tfp.open((pref + PEXT).data(), "rb");
        for(uint64_t hv, rank;tfp.read(hv) == sizeof(hv); pfp.write(rank)) {
            try {
                rank = map->operator[](hv).rank_;
#if MAKE_HISTOGRAM
                ++ranks[rank];
#endif
                // std::fprintf(stderr, "rank %" PRIu64 " for item %" PRIu64 " is item %zu\n", rank, hv, ++ranknum);
            } catch(const std::out_of_range &ex) {
                std::fprintf(stderr, "missing item %" PRIu64 " is number %zu\n", hv, ++ranknum);
                map->for_each([](const auto &x, const auto &y) {
                    LOG_DEBUG("%lu is the hash for %s\n", x, y.s_);
                });
                throw;
            }
        }
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
    FType ifp_;
    const char *prefix_;
    std::vector<char> buf_;
    std::vector<char> lastcharsvec_;
    std::vector<uint64_t> parsevec_;
    std::vector<uint64_t> *savec_;
    size_t             bi_;
    // Prime 1?
public:
    size_t skip, parse, words;
#if !NDEBUG
    size_t updates;
#endif

    khmap *map_;


    template<typename C>
    auto hash(const C &c) const {return h_.hash(c);}

    static constexpr uint64_t LARGE_PRIME = (1ull << 61);
    static constexpr uint64_t SMALL_PRIME = 109829;
    static constexpr uint64_t MEDIUMPRIME = 1999999973;
    static constexpr uint64_t WINDOW_MOD  = 100;

    HashPass(const HashPass &) = delete;
    HashPass(HashPass &&o) : h_(std::move(o.h_)), cstr_(std::move(o.cstr_)), q_(std::move(o.q_)),
        ifp_(std::move(o.ifp_)),
        prefix_(o.prefix_), buf_(std::move(o.buf_)),
        lastcharsvec_(std::move(o.lastcharsvec_)),
        parsevec_(std::move(o.parsevec_)),
        savec_(nullptr), 
        bi_(o.bi_),
        skip(o.skip), parse(o.parse), words(o.words)
    {
        std::swap(savec_, o.savec_);
        assert(o.savec_ == nullptr);
#if !NDEBUG
        updates = 0;
#endif
        make_map();
    }

    void make_map() {
        if(!map_) map_ = new khmap();
    }

    std::string str() const {
        std::string cs(cstr_.begin(), cstr_.end());
        std::string qs(q_.begin(), q_.end());
        char buf[2048];
        std::sprintf(buf, "HashPass:{[strings: %s|%s][fps: %s|%s][skip|parse|words:%zu|%zu|%zu]}",
                     cs.data(), qs.data(), ifp_.path(),
                     skip, parse, words);
        return buf;
    }

    ~HashPass() {
        if(map_) delete map_;
        if(savec_) delete savec_;
    }

    HashPass &operator|=(HashPass &&hp) {
        throw "a party";
        HashPass tmp(std::move(hp)); // To drop
        return *this;
    }

    HashPass(unsigned wsz, size_t nchunks=1, size_t chunknum=0, bool makesa=true, const char *pref="default_prefix", int compression=6):
        h_(wsz), i_(0), prefix_(pref), buf_(1<<14),
        savec_(makesa ? new std::vector<uint64_t>: nullptr),
        bi_(buf_.size()),
        map_(nullptr)
    {
        LOG_DEBUG("Starting hashpass with prefix=%s\n", pref);
        std::string safn = pref; safn += '.', safn += std::to_string(chunknum);
        safn += ".sa";
        std::string mode = "wb";
        if constexpr(FType::is_gz()) {
            if(compression == 0) mode = "wT";
            else mode += std::to_string(compression);
        }
#if !NDEBUG
        updates = 0;
#endif
    }
    void open_ifp(const char *path, const char *mode="rb") {
        ifp_.open(path, mode);
    }
    int nextchar() {
        if(__builtin_expect(bi_ == buf_.size(), 0)) {
            size_t n = ifp_.read(buf_.data(), buf_.size());
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
        LOG_DEBUG("About to fill from ifp = %s\n", ifp_.path().data());
        if(start) {
            ifp_.seek(start);
            int c;
            for(;(c = nextchar()) != EOF && q_.size() < unsigned(w());q_.push_back(c)) {
                h_.eat(c);
                q_.push_back(c);
                if(++skip == nelem + start) {LOG_WARNING("sequence too short\n"); return;}
            }
            if(q_.size() != unsigned(w())) {LOG_WARNING("could not fill. Returning early\n"); return;}
            while(skip < unsigned(w()) || h_.hashvalue % WINDOW_MOD) {
                if((c = nextchar()) == EOF) {LOG_WARNING("sequence too short2\n"); return;}
                if(++skip == nelem + start) {LOG_WARNING("sequence too short3\n"); return;}
                h_.update(q_.front(), c);
                q_.pop_front();
                q_.push_back(c);
            }
            parse = w();
            skip -= w();
            assert(unsigned(w()) == q_.size());
            cstr_ = q_;
        } else cstr_.push_front(util::Dollar);
        // std::fprintf(stderr, "size of cstr: %zu\n", cstr_.size());
        uint64_t pos = start;
        if(pos) pos += skip + w();
#if 0
#define print_cstr() {\
        std::string tmp(cstr_.begin(), cstr_.end());\
        std::fprintf(stderr, "tmp: %s. tmp.size(): %zu\n", tmp.data() + (start == 0), tmp.size());\
        print_as_values(cstr_);\
    }
#else
#define print_cstr()
#endif
        for(int c; (c = nextchar()) != EOF;) {
            ++parse;
            if(unlikely(q_.size() < unsigned(w()))) {
                h_.eat(static_cast<unsigned char>(c));
                q_.push_back(c);
                cstr_.push_back(c);
                print_cstr();
            } else {
                assert(q_.size() == unsigned(w()));
                h_.update(static_cast<unsigned char>(q_.front()), static_cast<unsigned char>(c));
                q_.pop_front();
                q_.push_back(c);
                cstr_.push_back(c);
                if(h_.hashvalue % WINDOW_MOD == 0) {
                    ++words;
                    update(&pos);
                    if(nelem && skip + parse == nelem + w()) break;
                }
            }
        }
        cstr_.insert(cstr_.end(), w(), util::Dollar);
        update(&pos);
        ifp_.close();
    }
    void update(uint64_t *pos) {
        auto hv = h_.hash(cstr_);
        parsevec_.push_back(hv);
        map_->insert(hv, cstr_);
        lastcharsvec_.push_back(cstr_[cstr_.size() - 1 - w()]);
        if(*pos==0) *pos = cstr_.size()-1;
        else       *pos += cstr_.size() - w();
        if(savec_) savec_->push_back(*pos);
        cstr_.erase(cstr_.begin(), cstr_.begin() + (cstr_.size() - w()));
        LOG_DEBUG("HashPass has had %zu updates.\n", ++updates);
    }
    auto w() const {return h_.n;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modlp(T v) {return v & LARGE_PRIME;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modmp(T v) {return v & MEDIUMPRIME;}
    template<typename T, typename=std::enable_if_t<std::is_integral_v<T>>>
    static constexpr auto modsp(T v) {return v & SMALL_PRIME;}
    size_t memory_usage() const {
        std::fprintf(stderr, "Calling memory_usage\n");
        return sizeof(*this) + buf_.size() + lastcharsvec_.size() +
        parsevec_.size() * sizeof(uint64_t) +
        (savec_ ? sizeof(*savec_) + savec_->size() * sizeof(uint64_t): 0) + sizeof(map_) +
        (map_ ? map_->memory_usage() : size_t(0));
    }
};


} // namespace dbt

#endif /* #ifndef PREFIX_FREE_PARSE_H__ */
