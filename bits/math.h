#ifndef OP_MATH_H
#define OP_MATH_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>


namespace op {
    // Returns base to the power of exp.
    template<int exp, class T>
    constexpr T pow(T base);

    // Returns floor(log(x, base)) if x > 0, otherwise it returns -1.
    template<int base, class T>
    constexpr int ilog(T x);

    // Returns the integer part of the square root of x.
    uint32_t isqrt(uint64_t x);
}



// Implementation.
namespace op {
    namespace detail {
        // Does the exponentiation by squaring for op::pow.
        template<int exp, class T>
        inline constexpr T pow_impl(T base, uint64_t result = 1) {
            return
                exp ? (
                    exp & 1 ? (
                        pow_impl<(exp >> 1)>(base * base, base * result)
                    ) : (
                        pow_impl<(exp >> 1)>(base * base, result)
                    )
                ) : (
                    result
                )
            ;
        }


        // Checks if base can be squared without overflow the type.
        template<class T>
        inline constexpr bool can_square(T base) {
            return base <= std::numeric_limits<T>::max() / base;
        }


        // Get the largest exponent x such that x is a power of two and pow(base, x) doesn't
        // overflow the type of base.
        template<class T, T base>
        inline constexpr int max_pot_exp(int result = 1) {
            return 
                can_square(base) ? (
                    // Despite the fact that this can never overflow we still have to check for
                    // overflow or the compiler will complain.
                    max_pot_exp<T, can_square(base) ? base * base : 1>(result * 2)
                ) : (
                    result
                )
            ;
        }


        // op::ilog implementation using a binary search.
        template<int base, class T, int i = max_pot_exp<T, base>()>
        inline constexpr int ilog_helper(T n, int x = 0) {
            return 
                i ? (
                    ilog_helper<base, T, i / 2>(
                        n >= pow<i, T>(base) ? n / pow<i, T>(base) : n,
                        n >= pow<i, T>(base) ? x + i : x
                    )
                ) : (
                    x
                )
            ;
        }
    }


    template<int exp, class T>
    inline constexpr T pow(T base) {
        static_assert(exp < 64, "pow exponents >= 64 can only overflow");

        return exp < 0 ? 1 / detail::pow_impl<-exp>(base) : detail::pow_impl<exp>(base);
    }


    template<int base, class T>
    inline constexpr int ilog(T x) {
        static_assert(!(base <= 0), "op::ilog is not useful for base <= 0");
        static_assert(base != 1, "op::ilog is not useful for base == 1");
        static_assert(std::is_integral<T>::value, "op::ilog only works on integer types");

        return x > 0 ? detail::ilog_helper<base>(x) : -1;
    }


    inline uint32_t isqrt(uint64_t x) {
        // For small values using native functions is fastest.
        // Where "small" means maximum integer accurate after float conversion.
        if (x <= 16785406ull) {
            return uint32_t(std::sqrt(float(x)));
        } if (x <= 4503599761588223ull) {
            return uint32_t(std::sqrt(double(x)));
        }

        // Newton's method, adapted from Hacker's Delight.
        // Since we know x > 4503599761588223 we can skip some branches.
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
