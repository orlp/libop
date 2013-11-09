#ifndef OP_RNG_H
#define OP_RNG_H

#include <cstdint>
#include <type_traits>

namespace op {
    class RNG {
        // Implementation of the Tyche-i RNG by Samuel Neves and Filipe Araujo.
        // Period is expected to be 2^127, with a very good uniform distribution.

    public:
        // If no seed is passed to the constructor it will use a platform-specific
        // seeding method, otherwise it will use the given seed. The index argument
        // is optional and safe to ignore, for more information on the use of it see
        // the Tyche-i paper (http://eden.dei.uc.pt/~sneves/pubs/2011-snfa2.pdf).
        RNG();
        RNG(uint64_t seed, uint32_t index = 0);

        // randuN returns a random integer between 0 and 2^N - 1, inclusive
        uint64_t randu64(); 
        uint32_t randu32();
        uint16_t randu16();
        uint8_t randu8();

        // randiN returns a random integer between -2^(N-1) and 2^(N-1), inclusive
        int64_t randi64();
        int32_t randi32();
        int16_t randi16();
        int8_t randi8();

        // fill buf with buflen random bytes
        void randbytes(unsigned char *buf, int buflen);

        // randrange returns a random integer between min and max, inclusive
        template <class T> T randrange(T min, T max);

        // shuffles the sequence given by (first, last) randomly
        template<class RandomIt> void shuffle(RandomIt first, RandomIt last);

    private:
        void init(uint64_t seed, uint32_t index);
        void mix();

        uint32_t a, b, c, d;
    };
}



// implementation

#include <algorithm>
#include <cstring>
#include <ctime>
#include <iterator>

#ifdef _WIN32
    #include <windows.h>
    #include <wincrypt.h>
#else
    #include <fstream>
#endif

// seeding and internal RNG
inline op::RNG::RNG() {
    static uint64_t default_seed = 0;
    static uint32_t default_seed_ctr = 0;

    while (default_seed == 0) {
        #ifdef _WIN32
            HCRYPTPROV prov;

            if (!CryptAcquireContext(&prov, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT))  {
                default_seed = std::time(0);
                break;
            }

            if (!CryptGenRandom(prov, 8, reinterpret_cast<unsigned char*>(&default_seed)))  {
                default_seed = std::time(0);
            }

            CryptReleaseContext(prov, 0);
        #else
            std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);

            if (!urandom.read(reinterpret_cast<unsigned char*>(&default_seed), 8)) {
                default_seed = std::time(0);
            }
        #endif
    }

    init(default_seed, default_seed_ctr);

    // if we exhausted our default seed, generate a new one next time
    if (++default_seed_ctr == 0) {
        default_seed = 0;
    }
}

inline op::RNG::RNG(uint64_t seed, uint32_t index) {
    init(seed, index);
}

inline void op::RNG::init(uint64_t seed, uint32_t index) {
    a = seed;
    b = seed >> 32;
    c = 2654435769u;
    d = 1367130551u ^ index;

    for (int i = 0; i < 20; ++i) {
        mix();
    }
}

inline void op::RNG::mix() {
    #define ROTL32(x, n) ((x << n) | (x >> (32 - n)))
    b = ROTL32(b, 25) ^ c;
    d = ROTL32(d, 24) ^ a;
    c -= d;
    a -= b;
    b = ROTL32(b, 20) ^ c;
    d = ROTL32(d, 16) ^ a;
    c -= d;
    a -= b;
    #undef ROTL32
}

// boring interfaces to the internal RNG
inline uint64_t op::RNG::randu64() {
    uint64_t tmp;
    mix(); tmp = b; mix();
    return (tmp << 32) | b;
}

inline uint32_t op::RNG::randu32() {
    mix();
    return b;
}

inline uint16_t op::RNG::randu16() {
    return randu32();
}

inline uint8_t op::RNG::randu8() {
    return randu32();
}

inline int64_t op::RNG::randi64() {
    uint64_t tmp = randu64();
    return reinterpret_cast<int64_t&>(tmp);
}

inline int32_t op::RNG::randi32() {
    uint32_t tmp = randu32();
    return reinterpret_cast<int32_t&>(tmp);
}

inline int16_t op::RNG::randi16() {
    uint16_t tmp = randu16();
    return reinterpret_cast<int16_t&>(tmp);
}

inline int8_t op::RNG::randi8() {
    uint8_t tmp = randu8();
    return reinterpret_cast<int8_t&>(tmp);
}

inline void op::RNG::randbytes(unsigned char *buf, int buflen) {
    if (!buf) return;

    while (buflen > 0) {
        mix();
        std::memcpy(buf, reinterpret_cast<unsigned char*>(&b), std::min(buflen, 4));
        buflen -= 4;
        buf += 4;
    }
}

template <class T>
inline T op::RNG::randrange(T min, T max) {
    if (min > max) {
        std::swap(min, max);
    }

    // find the size of the range and the lowest power of 2 that is greater or
    // equal to the given range (we actually find that power minus one)
    typename std::make_unsigned<T>::type range_size = max - min;
    typename std::make_unsigned<T>::type power_mask = range_size - 1;

    for (int i = 1; i <= (4 * sizeof(T)); i <<= 1) {
        power_mask |= power_mask >> i;
    }

    while (true) {
        typename std::make_unsigned<T>::type rand;

        randbytes(reinterpret_cast<unsigned char*>(&rand), sizeof(T));

        rand &= power_mask;

        if (rand <= range_size) {
            return min + rand;
        }
    }
}

template<class RandomIt>
void op::RNG::shuffle(RandomIt first, RandomIt last) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;

    // Fisher-Yates shuffle
    for (diff_t i = last - first - 1; i > 0; --i) {
        std::swap(first[i], first[randrange(diff_t(0), i)]);
    }
}

#endif