#ifndef PREFIX_FREE_PARSE_H__
#define PREFIX_FREE_PARSE_H__
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <deque>
#include <functional>
#include <cmath>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <vector>
#include "util.h"
#include "klib/khash.h"
#include "logutil.h"
#include "include/krw.h"
#include "clhash/include/clhash.h"
template<typename T> class TD;

#ifndef INLINE
#  if __GNUC__ || __clang__ || INLINE_OVERRIDE
#    define INLINE __attribute__((always_inline))
#  else
#    define INLINE
#  endif
#endif
#define C2B(x) (x ? "true": "false")

namespace dbt {
using u64 = uint64_t;
using u32 = uint32_t;
using ustring = std::basic_string<unsigned char>;
using namespace std::literals;

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

class khmap: khash_t(m) {
public:
    using value_type     = hit_t;
    using reference_type = value_type &;
    using pointer_type   = value_type *;
    using key_type       = u64;

    khmap() {std::memset(this, 0, sizeof(*this));}
    void free_map() {
        std::free(this->keys);
        std::free(this->flags);
        std::free(this->vals);
        std::memset(this, 0, sizeof(*this));
    }
    void free() {
        free_values();
        free_map();
    }
    auto get(u64 v) const {
        return kh_get(m, this, v);
    }
    auto       &val(u64 ki)       {return kh_val(this, ki);}
    const auto &val(u64 ki) const {return kh_val(this, ki);}
    key_type key(u64 ki) const {return kh_key(this, ki);}
    hit_t &operator[](u64 v) {
#define HIT_CORE\
        if(auto ki = get(v); unlikely(is_end(v)))\
            throw std::out_of_range(std::string("Missing key ") + std::to_string(v));\
        else\
            return val(ki);
        HIT_CORE
    }
    const hit_t &operator[](u64 v) const {
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
    size_t size() const {return kh_size((khash_t(m) *)this);}
    size_t capacity() const {return this->n_buckets;}
    bool is_end(khint_t ki) const {return ki == capacity();}
    bool exist(khint_t ki) const {return kh_exist(this, ki);}
    void print_all() const {
        for_each([](auto x, const hit_t &y) {
            std::fprintf(stderr, "String %s has hash %" PRIu64 "\n", y.s_, x);
        });
    }
    size_t memory_usage(bool add_strlens=false) const {
        size_t ret = sizeof(*this) + size() * (sizeof(value_type) + sizeof(key_type)) +
        __ac_fsize(this->n_buckets) * sizeof(uint32_t);
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
    auto put(u64 k) {
        int khr;
        return kh_put(m, this, k, &khr);
    }
    bool contains_key(u64 k) const {return !is_end(get(k));}
    // TODO: consider compressing strings by storing in 4-bits per character
    void insert(u64 v, const ustring &s) {
        insert(v, (const char *)s.data(), s.size());
    }
    void insert(u64 v, const char *s, size_t nelem) {
        if(auto ki = get(v); !is_end(ki)) {
            assert(this->vals[ki].s_);
            if(unlikely(std::strcmp(s, this->vals[ki].s_))) {
                char buf[4096];
                std::sprintf(buf, "Hash collision. Hash value: %" PRIu64 ". strings: '%s/%zu', '%s/%zu\n", v, s, std::strlen(s), this->vals[ki].s_, std::strlen(this->vals[ki].s_));
                throw std::runtime_error(buf);
            }
            if(unlikely(++this->vals[ki].occ_ == 0)) // This line both increments and checks! Do not remove.
                throw std::runtime_error("Overflow in occurrence count");
        } else {
            ki = put(v);
            char *s2 = static_cast<char *>(std::malloc(nelem + 1));
            if(!s2) throw std::bad_alloc();
            std::memcpy(s2, s, nelem);
            s2[nelem] = '\0';
            val(ki) = hit_t{const_cast<const char *>(s2), 1};
            std::set<char> set(s2, s2 + nelem);
#if 0
            assert(std::accumulate(s2, s2 + nelem, true, [](bool v, const signed char c) -> bool {
                return v && ((c >= 0) || c == util::Dollar);
            }) ||
                  !std::fprintf(stderr, "Contents of set: %s\n", stringify(set).data()));
#endif
        }
    }
    const khash_t(m) *map() const {return this;}
    void insert(u64 v, const char *s) {insert(v, s, std::strlen(s));}
    template<typename T>
    void insert(u64 v, const T &s) {
        std::string t(s.begin(), s.end());
        insert(v, t.data(), t.size());
    }
    khmap(khmap &&m) {std::memset(this, 0, sizeof(*this)); *this = std::move(m);}
    operator       khash_t(m) *()        {return reinterpret_cast<khash_t(m) *>(this);}
    operator const khash_t(m) *()  const {return reinterpret_cast<const khash_t(m) *>(this);}
        
    khmap &operator=(khmap &&m) {
        this->free();
        std::memcpy(this, &m, sizeof(m));
        std::memset(&m, 0, sizeof(m));
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
        for(khiter_t ki = 0; ki < kh_end(this); ++ki)
            if(exist(ki))
                func(key(ki), val(ki));
    }
    template<typename Func>
    void for_each(const Func &func) {
        for(khiter_t ki = 0; ki < kh_end(this); ++ki)
            if(exist(ki))
                func(key(ki), val(ki));
    }
    khmap &operator|=(khmap &&o) {
        assert_nonnull();
        o.assert_nonnull();
        o.for_each([&](const auto key, hit_t &h) {
            if(auto ki = this->get(key); ki == this->capacity()) {
                ki = this->put(key);
                val(ki) = h;
            } else {
                val(ki).occ_ += h.occ_;
                std::free(const_cast<char *>(h.s_));
            }
            std::memset(&h, 0, sizeof(h));
        });
        kh_clear(m, &o);
        return *this;
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


template<typename, typename> class HashPasser;


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

namespace constants {
    static constexpr u64 LARGE_PRIME = (1ull << 61);
    static constexpr u64 SMALL_PRIME = 109829;
    static constexpr u64 MEDIUMPRIME = 1999999973;
    static constexpr u64 WINDOW_MOD  = 100;
}

struct ResultSet {
    std::vector<char> lastcs_; // This is used to differentiate phrases which are also in $S$
    std::vector<u64>  parses_; // This contains hashes, which serve as identifiers for phrases, until they are replaced by their lexicographic rank in the sorted dictionary.
    std::vector<u64> sa_;
    khmap                    map_;
    size_t                 words_;
    size_t                parsed_;
    size_t               skipped_;
    unsigned                  id_; // Integral value reflecting which result subsection it is.
    u64 start_, stop_;
    bool make_sa_;
    ResultSet(u32 id, bool savec=false, u64 start=0, u64 stop=0): words_(0), parsed_(0), skipped_(0), id_(id), start_(start), stop_(stop), make_sa_(savec) {
        reserve_all(stop_ - start_);
        std::fprintf(stderr, "Made resultset with savec = %s\n", make_sa_ ? "true": "false");
    }
    ResultSet(ResultSet &&) = default;
    ResultSet(bool makesa=false) {
        std::memset(this, 0, sizeof(*this)); make_sa_ = makesa;
        std::fprintf(stderr, "Made resultset with boring constructor and = %s\n", make_sa_ ? "true": "false");
    }
    ResultSet &operator=(ResultSet &&) = default;
    void reserve_all(size_t size) {
        lastcs_.reserve(size); parses_.reserve(size); if(make_sa_) sa_.reserve(size);
    }
    void make_empty() {
#define DEL(x) decltype(x)().swap(x);
        DEL(lastcs_);
        DEL(parses_);
        DEL(sa_);
        kh_destroy(m, map_);
        // words_ = parsed_ = skipped_ = 0;
    }
    ResultSet &merge_and_destroy(ResultSet &&o) {
        if(o.id_ < id_) throw std::runtime_error("Attempting to merge results out of order. Abort!\n");
        if(!make_sa_ != !o.make_sa_) {
            char buf[2048];
            std::sprintf(buf, "ResultSets being merged where presence or absence of suffix array sampling do not agree (%s, %s)\n", C2B(make_sa_), C2B(o.make_sa_));
        }
        lastcs_.insert(lastcs_.end(), o.lastcs_.begin(), o.lastcs_.end());
        parses_.insert(parses_.end(), o.parses_.begin(), o.parses_.end());
        if(make_sa_) sa_.insert(sa_.end(), o.sa_.begin(), o.sa_.end());
        parsed_ += o.parsed_;
        words_ += o.words_;
        skipped_ += o.skipped_;

        if(map_.size() > o.map_.size()) // Merge into the larger map.
            std::swap_ranges((uint8_t *)&map_, (uint8_t *)&o.map_,  (uint8_t *)&o.map_ + sizeof(o.map_));
        map_ |= std::move(o.map_);
        return *this;
    }
};


template<typename PType>
static inline INLINE int nextchar(util::FpWrapper<PType> &ifp, size_t &bi, std::vector<char> &buf) {
    if(__builtin_expect(bi == buf.size(), 0)) {
        int n = ifp.bulk_read(buf.data(), buf.size()); /* To omit double buffering. */
        if(n <= 0) return EOF;
        if(n != static_cast<ssize_t>(buf.size())) buf.resize(n); // We're at the end,  so just shrink the buffer so that we don't eat garbage.
        bi = 0;
    }
    return buf[bi++];
}

namespace detail {
template<typename Hasher>
void update(ResultSet &rs, u64 *pos, ustring &cstr, const Hasher &h, u32 wsz) {
    auto hv = h(*(std::string *)&cstr);
    rs.parses_.push_back(hv);
    rs.map_.insert(hv, cstr);
    rs.lastcs_.push_back(cstr[cstr.size() - 1 - wsz]);
    if(*pos==0) *pos = cstr.size()-1;
    else       *pos += cstr.size() - wsz;
    if(rs.make_sa_) rs.sa_.push_back(*pos);
    cstr.erase(cstr.begin(), cstr.begin() + (cstr.size() - wsz));
}
}

/*
 Set 'nelem' to 0 to process a full file.
 */
template<typename Hasher=clhasher, typename PointerType=std::FILE *>
void perform_subwork(ResultSet &rs, const Hasher &hasher, u32 wsz, const char *path) {
    ustring cstr;
    krw::KRWindow krw(wsz);
    std::vector<char> buf(1 << 14);
    size_t bi = 0;
    util::FpWrapper<PointerType> ifp;
    ifp.open(path, "rb");
    rs.skipped_ = rs.parsed_ = rs.words_ = 0;
    u64 start = rs.start_, nelem = rs.stop_ - rs.start_;
    LOG_INFO("About to fill from ifp = %s, with [%zu-%zu)\n", path, start, start + nelem);
    if(start) {
        ifp.seek(start);
        int c;
        while((c = nextchar(ifp, bi, buf)) != EOF) {
            if(++rs.skipped_ == nelem + wsz) {LOG_WARNING("sequence too short\n"); return;}
            krw.eat(c);
            cstr.push_back(c);
            if(krw.hash % constants::WINDOW_MOD == 0 && rs.skipped_ >= wsz)
                break;
        }
#if 0
        if(krw.tot_char != wsz) {LOG_WARNING("could not fill. Returning early. fpos: %d. is eof: %s. total parsed: %d. k_: %u.\n", int(ifp.tell()), C2B(ifp.eof()), krw.tot_char, krw.k_); return;}
        while(rs.skipped_ < wsz || krw.hash % constants::WINDOW_MOD) {
            if((c = nextchar(ifp, bi, buf)) == EOF) {LOG_WARNING("sequence too short2\n"); return;}
            if(++rs.skipped_ == nelem + start) {LOG_WARNING("sequence too short3\n"); return;}
            krw.eat(c);
        }
#endif
        rs.parsed_   = wsz;
        rs.skipped_ -= wsz;
        if(cstr.size() > wsz) cstr.erase(0, cstr.size() - wsz);
    } else cstr.push_back(util::Dollar);
    u64 pos = start ? start + rs.skipped_ + wsz: 0;
    using detail::update;
    for(int c; likely((c = nextchar(ifp, bi, buf)) != EOF);) {
        ++rs.parsed_;
        if(cstr.size() >= wsz) {
            krw.eat(c);
            cstr.push_back(c);
            if(krw.hash % constants::WINDOW_MOD == 0) {
                ++rs.words_;
                update(rs, &pos, cstr, hasher, wsz);
                if(nelem && rs.skipped_ + rs.parsed_ == nelem + wsz) {
                    LOG_INFO("Number of expected elements finished.\n");
                }
            }
        } else krw.eat(static_cast<unsigned char>(c)), cstr.push_back(c);
    }
    cstr.insert(cstr.end(), wsz, util::Dollar);
    update(rs, &pos, cstr, hasher, wsz);
}

template<typename PointerType=std::FILE *, typename Hasher=std::hash<std::string>>
class HashPasser {
    using FType = util::FpWrapper<PointerType>;
    u32 wsz_;
    Hasher h_;
    std::vector<ResultSet> results_;
    std::string path_;
public:


    HashPasser(const char *path, unsigned wsz=10, int nthreads=0, bool makesa=true):
        wsz_(wsz), path_(path)
    {
        if(nthreads < 1)
            nthreads = std::thread::hardware_concurrency(); // I want all cores and threads you have.
        make_subs(path, nthreads, makesa);
    }
    void make_subs(const char *path, unsigned nthreads, bool makesa) {
        const size_t fsz = util::get_fsz<PointerType>(path), per_chunk = std::ceil(double(fsz) / nthreads);
        size_t start = 0, stop;
        results_.reserve(nthreads);
        while(results_.size() < unsigned(nthreads)) {
            stop = std::min(fsz, start + per_chunk);
            results_.emplace_back(results_.size(), makesa, start, stop);
            start = stop;
        }
    }
    void seed(uint64_t newseed) {h_.seed(newseed);} //
    void run() {
        // #pragma omp parallel
        for(size_t i = 0; i < results_.size(); ++i) {
            perform_subwork(results_[i], h_, wsz_, path_.data());
        }
        // TODO: improve the 'reduce' portion of this.
        // It really should be done in parallel in a divide and conquer tree of logarithmic depth.
        std::fprintf(stderr, "About to merge and destroy from other chunks %zu\n", results_.size());
        for(size_t i = 1; i < results_.size(); ++i)
            results_[0].merge_and_destroy(std::move(results_[i]));
        results_.erase(results_.begin() + 1, results_.end());
        std::fprintf(stderr, "Mrged and destroyed from other chunks %zu\n", results_.size());
    }
    HashPasser(HashPasser &o) = default;
    HashPasser(const HashPasser &) = delete;

    ~HashPasser() {
    }

    auto w() const {return h_.n;}
};

#undef C2B

} // namespace dbt

#endif /* #ifndef PREFIX_FREE_PARSE_H__ */
