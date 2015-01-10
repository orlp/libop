#ifndef OP_MATH_H
#define OP_MATH_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <set>
#include <iterator>

#include "random.h"


namespace op {
    // Returns base to the power of exp.
    template<int exp, class T>
    constexpr T pow(T base);

    // Returns floor(log(x, base)) if x > 0, otherwise it returns -1.
    template<int base, class T>
    constexpr int ilog(T x);

    // Returns the integer part of the square root of x.
    template<class T>
    inline uint32_t isqrt(T x);

    // Returns x < y, doing it correctly even if the signedness differs between T and U.
    template<class T, class U>
    constexpr bool safe_less(T x, U y);

    // Returns x > y, doing it correctly even if the signedness differs between T and U.
    template<class T, class U>
    constexpr bool safe_greater(T x, U y);

    // Returns x <= y, doing it correctly even if the signedness differs between T and U.
    template<class T, class U>
    constexpr bool safe_less_equal(T x, U y);

    // Returns x >= y, doing it correctly even if the signedness differs between T and U.
    template<class T, class U>
    constexpr bool safe_greater_equal(T x, U y);
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


    template<class T>
    inline uint32_t isqrt(T x) {
        static_assert(std::is_integral<T>::value, "op::isqrt only works on integer types");

        if (x < 0) throw std::domain_error("op::isqrt");
        return op::isqrt<uint64_t>(x);
    }

    template<>
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

    
    template<class T, class U>
    inline uint64_t gcd(T a_, U b_) {
        static_assert(std::is_integral<T>::value && std::is_integral<U>::value,
                      "op::gcd only works on integer types");

        uint64_t a = a_ >= 0 ? a_ : -a_;
        uint64_t b = b_ >= 0 ? b_ : -b_;
        if (a == b) return a;
        
        while (b > 0) {
            uint64_t tmp = a;
            a = b;
            b = tmp % b;
        }

        return a;
    }


    template<class T, class U>
    inline uint64_t lcm(T a_, U b_) {
        static_assert(std::is_integral<T>::value && std::is_integral<U>::value,
                      "op::lcm only works on integer types");

        uint64_t a = a_ >= 0 ? a_ : -a_;
        uint64_t b = b_ >= 0 ? b_ : -b_;

        return (a / gcd(a, b)) * b;
    }


    template<class OutIter>
    inline void primesbelow(uint64_t limit, OutIter out) {
        // http://stackoverflow.com/questions/4643647/fast-prime-factorization-module
        if (limit > 2) *out++ = 2;
        if (limit > 3) *out++ = 3;
        if (limit <= 5) return;

        uint64_t correction = limit % 6 > 1;
        uint64_t wheels[6] = { limit, limit - 1, limit + 4, limit + 3, limit + 2, limit + 1 };
        uint64_t n = wheels[limit % 6];

        std::vector<bool> sieve(n / 3, true);
        sieve[0] = false;

        for (uint64_t i = 0, upper = op::isqrt(n)/3; i <= upper; ++i) {
            if (sieve[i]) {
                uint64_t k = (3*i + 1) | 1;
                for (uint64_t j = (k*k) / 3;                   j < n/3; j += 2*k) sieve[j] = false;
                for (uint64_t j = (k*k + 4*k - 2*k*(i & 1))/3; j < n/3; j += 2*k) sieve[j] = false;
            }
        }

        for (uint64_t i = 1; i < n / 3 - correction; ++i) {
            if (sieve[i]) *out++ = (3 * i + 1) | 1;
        }
    }


    namespace detail {
        // Computes (a + b) % m, assumes a < m, b < m.
        uint64_t addmod64(uint64_t a, uint64_t b, uint64_t m) {
            if (b >= m - a) return a - m + b;
            return a + b;
        }

        // Computes (a*b) % m safely, considering overflow. Requires b < m;
        uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
            // No overflow possible.
            if (b <= std::numeric_limits<uint64_t>::max() / a) return (a*b) % m;

            uint64_t res = 0;
            while (a != 0) {
                if (a & 1) res = addmod64(res, b, m);
                a >>= 1;
                b = addmod64(b, b, m);
            }

            return res;
        }
    }


    uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
        uint64_t r = 1;

        b %= m;
        while (e) {
            if (e % 2 == 1) r = detail::mulmod64(r, b, m);
            e >>= 1;
            b = detail::mulmod64(b, b, m);
        }

        return r;
    }


    template<class T>
    inline bool isprime(T x) {
        static_assert(std::is_integral<T>::value, "op::isprime only works on integer types");

        if (x < 2) return false;
        return op::isprime<uint64_t>(x);
    }


    template<>
    inline bool isprime(uint64_t n) {
        constexpr int max_smallprimeset = 100000;
        static std::set<int> smallprimeset;

        if (smallprimeset.size() == 0) {
            op::primesbelow(max_smallprimeset, std::inserter(smallprimeset, smallprimeset.begin()));
        }

        if (n <= 3) return n >= 2;
        if (n % 2 == 0) return false;
        if (n < max_smallprimeset) return smallprimeset.count(n);

        uint64_t d = n - 1;
        int s = -1;
        while (d % 2 == 0) {
            d /= 2;
            s += 1;
        }

        uint64_t bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
        for (auto base : bases) {
            uint64_t x = op::powmod(base, d, n);

            if (x == 1 || x == n - 1) continue;
            for (int i = 0; i < s; ++i) {
                x = detail::mulmod64(x, x, n);

                if (x == 1) return false;
                if (x == n - 1) goto next_base;
            }

            return false;
        next_base: ;
        }

        return true;
    }
    

    namespace detail {
        uint64_t pollard_brent(uint64_t n) {
            static std::mt19937_64 rng((op::random_device()()));

            uint64_t y = op::randint<uint64_t>(1, n-1, rng);
            uint64_t c = op::randint<uint64_t>(1, n-1, rng);
            uint64_t m = op::randint<uint64_t>(1, n-1, rng);

            uint64_t g, r, q, x, ys;
            g = r = q = 1;

            while (g == 1) {
                x = y;

                for (uint64_t i = 0; i < r; ++i) {
                    y = detail::addmod64(detail::mulmod64(y, y, n), c, n);
                }

                for (uint64_t k = 0; k < r && g == 1; k += m) {
                    ys = y;

                    for (uint64_t i = 0; i < std::min(m, r-k); ++i) {
                        y = detail::addmod64(detail::mulmod64(y, y, n), c, n);
                        q = detail::mulmod64(q, x < y ? y-x : x-y, n);
                    }

                    g = op::gcd(q, n);
                }

                r *= 2;
            }

            if (g == n) {
                do {
                    ys = detail::addmod64(detail::mulmod64(ys, ys, n), c, n);
                    g = op::gcd(x < ys ? ys-x : x-ys, n);
                } while (g == 1);
            }

            return g;
        }
    }


    inline std::vector<uint64_t> primefactors(uint64_t n) {
        static std::vector<int> smallprimes;
        if (smallprimes.size() == 0) op::primesbelow(1000, std::back_inserter(smallprimes));

        std::vector<uint64_t> factors;
        for (uint64_t checker : smallprimes) {
            while (n % checker == 0) {
                factors.push_back(checker);
                n /= checker;
            }

            if (checker > n) break;
        }

        std::vector<uint64_t> to_factor = {n};
        while (to_factor.size()) {
            n = to_factor.back();
            to_factor.pop_back();

            if (n == 1) continue;
            if (op::isprime(n)) {
                factors.push_back(n);
                continue;
            }

            // Get (not necessarily prime) factor of n.
            uint64_t factor = detail::pollard_brent(n);
            to_factor.push_back(factor);
            to_factor.push_back(n / factor);
        }

        return factors;
    }
    




    namespace detail {
        template<class T, class U, bool=std::is_signed<T>::value, bool=std::is_signed<U>::value>
        struct safe_less_helper;

        template<class T, class U, bool same_sign>
        struct safe_less_helper<T, U, same_sign, same_sign> {
            static constexpr bool eval(T x, U y) { return x < y; }
        };

        template<class T, class U>
        struct safe_less_helper<T, U, true, false> {
            static constexpr bool eval(T x, U y) {
                return x < 0 || static_cast<typename std::make_unsigned<T>::type>(x) < y;
            }
        };

        template<class T, class U>
        struct safe_less_helper<T, U, false, true> {
            static constexpr bool eval(T x, U y) {
                return y > 0 && x < static_cast<typename std::make_unsigned<U>::type>(y);
            }
        };
    }

    template<class T, class U>
    inline constexpr bool safe_less(T x, U y) {
        // TODO: allow float-integer comparison
        static_assert(std::is_integral<T>::value == std::is_integral<U>::value,
                      "op::safe_less only works if both types are integer or float");

        return detail::safe_less_helper<T, U>::eval(x, y);
    }

    template<class T, class U>
    inline constexpr bool safe_greater(T x, U y) { return safe_less(y, x); }

    template<class T, class U>
    inline constexpr bool safe_less_equal(T x, U y) { return !safe_less(y, x); }

    template<class T, class U>
    inline constexpr bool safe_greater_equal(T x, U y) { return !safe_less(x, y); }
}

#endif
