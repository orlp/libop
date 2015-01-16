#ifndef OP_INTRIN_H
#define OP_INTRIN_H

namespace op {
    // Multiplies two unsigned 64 bit numbers. Returns a pair of 64 bit numbers containing the
    // result (high, low).
    inline std::pair<uint64_t, uint64_t> mulu64(uint64_t a, uint64_t b);
    
    // Returns the remainder after dividing the unsigned 128 bit number (hi, lo) by m.
    inline uint64_t modu128_64(uint64_t hi, uint64_t lo, uint64_t m);

    // Returns (a*b) % m without risk of overflow.
    inline uint64_t mulmodu64(uint64_t a, uint64_t b, uint64_t m);
}



// Implementation.
namespace op {
    inline std::pair<uint64_t, uint64_t> mulu64(uint64_t a, uint64_t b) {
        #if defined(__GNUC__) && defined(__x86_64__)
            uint64_t h, l;
            asm("mulq %3"
                : "=a"(l),"=d"(h)
                : "a"(a), "rm"(b)
                : "cc");
            return std::make_pair(h, l);
        #else
            #error no mulu64 implementation
        #endif

        // TODO: fallback implementation
    }


    inline uint64_t mulmodu64(uint64_t a, uint64_t b, uint64_t m) {
        #if defined(__GNUC__) && defined(__x86_64__)
            uint64_t q, r;
            asm("mulq %3;"
                "divq %4;"
                : "=a"(q), "=d"(r)
                : "a"(a), "d"(b), "rm"(m)
                : "cc");
            return r;
        #else
            a %= m;
            b %= m;

            // No overflow possible.
            if (a == 0) return b;
            if (b <= std::numeric_limits<uint64_t>::max() / a) return (a*b) % m;

            uint64_t res = 0;
            while (a != 0) {
                if (a & 1) {
                    if (b >= m - res) res -= m;
                    res += b;
                }

                a >>= 1;
                if (b >= m - b) b += b - m;
                else            b += b;
            }

            return res;
        #endif
    }


    inline uint64_t modu128_64(uint64_t hi, uint64_t lo, uint64_t m) {
        #if defined(__GNUC__) && defined(__x86_64__)
            uint64_t q, r;
            asm("divq %4"
                : "=a"(q),"=d"(r)
                : "a"(lo), "d" (hi), "rm"(m)
                : "cc");
            return r;
        #else
            #error no modu128_64 implementation
        #endif

        // TODO: fallback implementation
    }
}

#endif
