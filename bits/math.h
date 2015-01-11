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
        inline uint64_t addmod64(uint64_t a, uint64_t b, uint64_t m) {
            if (b >= m - a) return a - m + b;
            return a + b;
        }

        // Computes (a || b) % m.
        inline uint64_t mod64(uint64_t a, uint64_t b, uint64_t m) {
            #if defined(__GNUC__) && defined(__x86_64__)
                uint64_t q, r;
                asm("divq %4"
                    : "=a"(q),"=d"(r)
                    : "a"(b), "d" (a), "rm"(m)
                    : "cc");
                return r;
            #else
                #error no mod64 implementation
            #endif

            // TODO: fallback implementation
        }


        inline std::pair<uint64_t, uint64_t> mul64(uint64_t a, uint64_t b) {
            #if defined(__GNUC__) && defined(__x86_64__)
                uint64_t h, l;
                asm("mulq %3"
                    : "=a"(l),"=d"(h)
                    : "a"(a), "rm"(b)
                    : "cc");
                return std::make_pair(h, l);
            #else
                #error no mul64 implementation
            #endif

            // TODO: fallback implementation
        }


        // Finds 2^-64 mod m and (-m)^-1 mod m for odd m (hacker's delight).
        inline std::pair<uint64_t, uint64_t> mont_modinv(uint64_t m) {
            uint64_t a = 1ull << 63;
            uint64_t u = 1;
            uint64_t v = 0;

            while (a > 0) {
                a = a >> 1;
                if ((u & 1) == 0) {
                    u = u >> 1; v = v >> 1;
                } else {
                    u = ((u ^ m) >> 1) + (u & m);
                    v = (v >> 1) + (1ull << 63);
                }
            }

            return std::make_pair(u, v);
        }


        // Computes aR * bR mod N with R = 2**64.
        inline uint64_t montmul64(uint64_t a, uint64_t b, uint64_t N, uint64_t Nneginv) {
            uint64_t Th, Tl, m, mNh, mNl, th;

            std::tie(Th, Tl) = mul64(a, b);
            m = Tl * Nneginv;
            std::tie(mNh, mNl) = mul64(m, N);

            bool lc = Tl + mNl < Tl;
            th = Th + mNh + lc;
            bool hc = (th < Th) || (th == Th && lc);

            if (hc > 0 || th >= N) th = th - N;

            return th;
        }


        // Computes (a*b) % m safely, considering overflow. Requires b < m;
        inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
            #if defined(__GNUC__) && defined(__x86_64__)
                uint64_t q, r;
                asm("mulq %3;"
                    "divq %4;"
                    : "=a"(q), "=d"(r)
                    : "a"(a), "d"(b), "rm"(m)
                    : "cc");
                return r;
            #else
                // No overflow possible.
                if (a == 0) return b;
                if (b <= std::numeric_limits<uint64_t>::max() / a) return (a*b) % m;

                uint64_t res = 0;
                while (a != 0) {
                    if (a & 1) res = addmod64(res, b, m);
                    a >>= 1;
                    b = addmod64(b, b, m);
                }

                return res;
            #endif
        }
    }


    inline uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
        uint64_t r = 1;

        b %= m;
        while (e) {
            if (e % 2 == 1) r = detail::mulmod64(r, b, m);
            b = detail::mulmod64(b, b, m);
            e >>= 1;
        }

        return r;
    }

    namespace detail {
        template<class Iter>
        inline bool miller_rabin(uint64_t n, Iter base_begin, Iter base_end) {
            // We only ever have to compare to 1 and n-1, instead of converting back and forth,
            // compute once and compare in Montgomery form.
            uint64_t nneginv = detail::mont_modinv(n).second;
            uint64_t mont1 = detail::mod64(1, 0, n);
            uint64_t montn1 = detail::mod64(n-1, 0, n);

            uint64_t d = n - 1;
            int s = -1;
            while (d % 2 == 0) {
                d /= 2;
                s += 1;
            }

            for (Iter it = base_begin; it != base_end; ++it) {
                uint64_t x = mont1;
                uint64_t b = detail::mod64(*it % n, 0, n);

                uint64_t e = d;
                while (e) {
                    if (e % 2 == 1) x = detail::montmul64(x, b, n, nneginv);
                    b = detail::montmul64(b, b, n, nneginv);
                    e >>= 1;
                }

                if (x == mont1 || x == montn1) continue;

                for (int i = 0; i < s; ++i) {
                    x = detail::montmul64(x, x, n, nneginv);
                    if (x == mont1) return false;
                    if (x == montn1) goto next_base;
                }

                return false;
            next_base: ;
            }

            return true;
        }
    }


    template<class T>
    inline bool isprime(T x) {
        static_assert(std::is_integral<T>::value, "op::isprime only works on integer types");

        if (x < 2) return false;
        return op::isprime<uint64_t>(x);
    }


    template<>
    inline bool isprime(uint64_t n) {
        // Miller-rabin with a check for small primes first.
        constexpr int max_smallprimeset = 100000;
        static std::set<int> smallprimeset;

        if (smallprimeset.size() == 0) {
            op::primesbelow(max_smallprimeset, std::inserter(smallprimeset, smallprimeset.begin()));
        }

        if (n <= 3) return n >= 2;
        if (n % 2 == 0) return false;
        if (n < max_smallprimeset) return smallprimeset.count(n);

        uint64_t s1[] = {9345883071009581737ull};
        uint64_t s2[] = {336781006125ull, 9639812373923155ull};
        uint64_t s3[] = {4230279247111683200ull, 14694767155120705706ull, 16641139526367750375ull};
        uint64_t s4[] = {2ull, 141889084524735ull, 1199124725622454117ull, 11096072698276303650ull};
        uint64_t s5[] = {2ull, 4130806001517ull, 149795463772692060ull, 186635894390467037ull,
                         3967304179347715805ull};
        uint64_t s6[] = {2ull, 123635709730000ull, 9233062284813009ull, 43835965440333360ull,
                         761179012939631437ull, 1263739024124850375ull};
        uint64_t s7[] = {2ull, 325ull, 9375ull, 28178ull, 450775ull, 9780504ull, 1795265022ull};

        if (n < 341531) return detail::miller_rabin(n, std::begin(s1), std::end(s1));
        if (n < 1050535501) return detail::miller_rabin(n, std::begin(s2), std::end(s2));
        if (n < 350269456337) return detail::miller_rabin(n, std::begin(s3), std::end(s3));
        if (n < 55245642489451) return detail::miller_rabin(n, std::begin(s4), std::end(s4));
        if (n < 7999252175582851) return detail::miller_rabin(n, std::begin(s5), std::end(s5));
        if (n < 585226005592931977) return detail::miller_rabin(n, std::begin(s6), std::end(s6));
        return detail::miller_rabin(n, std::begin(s7), std::end(s7));
    }

    

    namespace detail {
        inline uint64_t pollard_brent(uint64_t n) {
            // Richard P. Brent (1980) – An Improved Monte Carlo Factorization Algorithm
            //
            // A cool optimization is used here. Firstly we use the much faster Montgomery reduction
            // for modular multiplication. However, note that:
            //
            //     a = b (mod n)   =>    gcd(a, n) = gcd(b, n)
            //
            // Since our result is always derived through a gcd involving n, and all logic and
            // branching is done on the result of such a gcd, it means we'll never have to convert
            // between Montgomery form and regular integers. In other words, Pollard-Brent just
            // works if you change all multiplications with Montgomery multiplication!
            //
            // The only precomputation needed is (-N)*2^64 (mod n)
            
            // Really simple LCG with a twist. Don't know how it works, but it works and is faster
            // than "better" RNGs for some reason.
            static uint64_t rng = 0xdeafbeef;
            uint64_t a = rng*6364136223846793005ull + 1442695040888963407ull;
            uint64_t b = a*6364136223846793005ull + 1442695040888963407ull;
            rng = (a+b) ^ (a*b);

            uint64_t y = 1 + a % (n - 1);
            uint64_t c = 1 + b % (n - 1);
            uint64_t m = 100;
            
            uint64_t nneginv = detail::mont_modinv(n).second;

            uint64_t g, r, q, x, ys;
            q = r = 1;

            do {
                x = y;
                for (uint64_t i = 0; i < r; ++i) {
                    y = detail::addmod64(detail::montmul64(y, y, n, nneginv), c, n);
                }

                for (uint64_t k = 0; k < r && g == 1; k += m) {
                    ys = y;
                    for (uint64_t i = 0; i < std::min(m, r-k); ++i) {
                        y = detail::addmod64(detail::montmul64(y, y, n, nneginv), c, n);
                        q = detail::montmul64(q, x < y ? y-x : x-y, n, nneginv);
                    }

                    g = op::gcd(q, n);
                }

                r *= 2;
            } while (g == 1);

            if (g == n) {
                do {
                    ys = detail::addmod64(detail::montmul64(ys, ys, n, nneginv), c, n);
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

            // Get a (not necessarily prime) factor of n.
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
