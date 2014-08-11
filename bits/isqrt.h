#ifndef OP_ISQRT_H
#define OP_ISQRT_H

#include <cstdint>
#include <cmath>

namespace op {
    /* returns the integer part of the square root of x */
    uint32_t isqrt(uint64_t x) {
        /* for small values using native functions is fastest */
        /* where "small" means maximum integer accurate after float conversion */
        if (x <= 16785406ull) {
            return uint32_t(sqrtf(float(x)));
        } if (x <= 4503599761588223ull) {
            return uint32_t(sqrt(double(x)));
        }

        /* Newton's method, adapted from Hacker's Delight */
        /* since we know x > 4503599761588223 we can skip some branches */
        uint64_t s, x1, g0, g1;
        x1 = (x - 1) >> 48;
        s = 25;
        if (x1 > 255) {s = s + 4; x1 = x1 >> 8;}
        if (x1 > 15)  {s = s + 2; x1 = x1 >> 4;}
        if (x1 > 3)   {s = s + 1;}

        g0 = 1 << s;
        g1 = (g0 + (x >> s)) >> 1;
        while (g1 < g0) {
            g0 = g1;
            g1 = (g0 + x/g0) >> 1;
        }

        return g0;
    }
}

#endif