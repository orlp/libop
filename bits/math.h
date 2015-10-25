#ifndef OP_MATH_H
#define OP_MATH_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "intrin.h"
#include "type_traits.h"


namespace op {
    // Returns base to the power of exp.
    template<int exp, class T>
    constexpr T pow(T base);
    
    // Returns base to the power of exp modulo mod.
    uint64_t powmod(uint64_t base, uint64_t exp, uint64_t mod);

    // Returns base to the power of exp.
    int64_t ipow(int32_t base, uint8_t exp);

    // Returns floor(log(x, base)) if x > 0, otherwise it returns -1.
    template<int base, class T>
    constexpr int ilog(T x);

    // Returns the integer part of the square root of integer x.
    template<class T>
    uint32_t isqrt(T x);

    // Hypothenuse in N dimensions.
    template<class... T>
    constexpr op::common_floating_point_type_t<T...> hypot(T... args);

    // Returns the greatest common divisor of a and b.
    template<class T, class U>
    uint64_t gcd(T a, U b);

    // Returns the least common multiple of a and b.
    template<class T, class U>
    uint64_t lcm(T a, U b);

    // Streams every prime less than limit as an uin64_t into out.
    template<class OutIter>
    void primes_below(uint64_t limit, OutIter out);

    // Returns true if n is prime.
    template<class T>
    bool is_prime(T n);

    // Returns a vector containing the prime factors of n.
    std::vector<uint64_t> prime_factors(uint64_t n);

    // Returns a map with the prime factors of n as the key, and how many times the prime factor
    // occurs in the factorization of n as the value.
    std::map<uint64_t, int> factorization(uint64_t n);    

    // Returns x < y, with correct results for all arithmetic types T and U, even if signedness
    // differs, or comparing floating point numbers to integers.
    template<class T, class U>
    constexpr bool safe_less(const T& x, const U& y);

    // Returns x > y, analogue to safe_less.
    template<class T, class U>
    constexpr bool safe_greater(const T& x, const U& y);

    // Returns x <= y, analogue to safe_less.
    template<class T, class U>
    constexpr bool safe_less_equal(const T& x, const U& y);

    // Returns x >= y, analogue to safe_less.
    template<class T, class U>
    constexpr bool safe_greater_equal(const T& x, const U& y);
}



// Implementation.
namespace op {
    namespace detail {
        // Does the exponentiation by squaring for op::pow.
        template<int exp, class T>
        inline constexpr T pow_impl(T base, uint64_t result = 1) {
            return
                exp ? (
                    (exp & 1) ? (
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


    inline int64_t ipow(int32_t base, uint8_t exp) {
        static const uint8_t highest_bit_set[] = {
              0,   1,   2,   2,   3,   3,   3,   3,
              4,   4,   4,   4,   4,   4,   4,   4,
              5,   5,   5,   5,   5,   5,   5,   5,
              5,   5,   5,   5,   5,   5,   5,   5,
              6,   6,   6,   6,   6,   6,   6,   6,
              6,   6,   6,   6,   6,   6,   6,   6,
              6,   6,   6,   6,   6,   6,   6,   6,
              6,   6,   6,   6,   6,   6,   6
        };

        if (exp >= 63) {
            if (base == 1) return 1;
            if (base == -1) return 1 - 2 * (exp & 1);

            throw std::domain_error("op::ipow");
        }

        uint64_t result = 1;
        switch (highest_bit_set[exp]) {
        case 6:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 5:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 4:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 3:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 2:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 1:
            if (exp & 1) result *= base;
        default:
            return result;
        }
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
        if (x <= 16785406ul) return uint32_t(std::sqrt(float(x)));
        return uint32_t(std::sqrt(double(x)));
    }

    template<>
    inline uint32_t isqrt(uint64_t x) {
        // For small values using native functions is fastest.
        // Where "small" means maximum integer accurate after float conversion.
        if (x <=         16785406ull) return uint32_t(std::sqrt(float(x)));
        if (x <= 4503599761588223ull) return uint32_t(std::sqrt(double(x)));

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

    namespace detail {
        template<class T>
        inline constexpr T hypot_impl_var(T s, T t) { return t * std::sqrt(s); }

        template<class T, class... A>
        inline constexpr T hypot_impl_var(T s, T t, T x, A... args) {
            return t < x  ? hypot_impl_var(T(1) + s * ((t/x)*(t/x)), x, args...) :
                   t != 0 ? hypot_impl_var(s + (x/t)*(x/t), t, args...) :
                            hypot_impl_var(s, t, args...);
        }

        template<class T>
        inline constexpr T hypot_impl(T x, T y) {
            return x < y  ? x * std::sqrt(T(1) + (x/y)*(x/y)) :
                   x != 0 ? x * std::sqrt(T(1) + (y/x)*(y/x)) :
                            0;
        }

        template<class T, class... A>
        inline constexpr T hypot_impl(T x, T y, A... args) {
            return x < y  ? hypot_impl_var(T(1) + (x/y)*(x/y), y, args...) :
                   x != 0 ? hypot_impl_var(T(1) + (y/x)*(y/x), x, args...) :
                            hypot_impl_var(T(1), T(0), args...);
        }
    }

    template<class... T>
    inline constexpr op::common_floating_point_type_t<T...> hypot(T... args) {
        return detail::hypot_impl(std::fabs(op::common_floating_point_type_t<T...>(args))...);
    }
    
    template<class T, class U>
    inline uint64_t gcd(T a_, U b_) {
        static_assert(std::is_integral<T>::value && std::is_integral<U>::value,
                      "op::gcd only works on integer types");

        uint64_t a = a_ >= 0 ? a_ : -a_;
        uint64_t b = b_ >= 0 ? b_ : -b_;
        if (!a || !b) return a | b;
        
        #if defined(__GNUC__) && defined(__x86_64__)
            int shift = __builtin_ctzll(a | b);
            a >>= __builtin_ctzll(a);
        
            while (b) {
                b >>= __builtin_ctzll(b);

                if (a < b) b -= a;
                else {
                    uint64_t t = a - b;
                    a = b;
                    b = t;
                }
            }

            a <<= shift;
        #else
            while (b > 0) {
                uint64_t t = a;
                a = b;
                b = t % b;
            }
        #endif

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
    inline void primes_below(uint64_t limit, OutIter out) {
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
        
        // Computes (a - b) % m, assumes a < m, b < m.
        inline uint64_t submod64(uint64_t a, uint64_t b, uint64_t m) {
            if (a < b) return m + a - b;
            return a - b;
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


        // Computes cR such that cR = aR * bR (mod n) with R = 2**64.
        inline uint64_t montmul64(uint64_t aR, uint64_t bR, uint64_t n, uint64_t nneginv) {
            uint64_t Th, Tl, m, mnh, mnl, th;

            std::tie(Th, Tl) = op::mulu64(aR, bR);
            m = Tl * nneginv;
            std::tie(mnh, mnl) = op::mulu64(m, n);

            bool lc = Tl + mnl < Tl;
            th = Th + mnh + lc;
            bool hc = (th < Th) || (th == Th && lc);

            if (hc > 0 || th >= n) th = th - n;

            return th;
        }
    }


    inline uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
        uint64_t r = 1;

        b %= m;
        while (e) {
            if (e % 2 == 1) r = op::mulmodu64(r, b, m);
            b = op::mulmodu64(b, b, m);
            e >>= 1;
        }

        return r;
    }

    namespace detail {
        // Assumes n is odd, does a strong probable prime test.
        inline bool strong_probable_prime(uint64_t n, uint64_t base,
                                          uint64_t mont1, uint64_t montn1, uint64_t nneginv) {
            uint64_t d = n - 1;
            int s = 0;
            while ((d & 1) == 0) {
                s += 1;
                d >>= 1;
            }

            // x = b^d (mod n)
            uint64_t x = mont1;
            uint64_t b = op::divu128_64(base % n, 0, n).second;
            uint64_t e = d;

            if (b == 0) return true;

            while (e) {
                if (e & 1) x = detail::montmul64(x, b, n, nneginv);
                b = detail::montmul64(b, b, n, nneginv);
                e >>= 1;
            }

            if (x == mont1 || x == montn1) return true;


            for (int r = 1; r < s; ++r) {
                x = detail::montmul64(x, x, n, nneginv);
                if (x == mont1) return false; // Early exit, 1^2 = 1 != n - 1.
                if (x == montn1) return true;
            }

            return false;
        }

        // Returns the Jacobi symbol of a over n. Assumes n odd.
        inline int jacobi(int64_t a_, uint64_t n) {
            int r = 1;

            // Handle negative a.
            uint64_t a;
            if (a_ < 0) {
                if (n % 4 == 3) r = -r;
                a = -a_;
            } else a = a_;

            while (true) {
                a %= n;

                // Divisibility by 4 can only cause double sign flips = no sign flips.
                while (a && a % 4 == 0) a /= 4; 
                if (a % 2 == 0) {
                    a /= 2;
                    if ((n % 8) == 3 || (n % 8) == 5) r = -r;
                }
                
                if (a == 0) return 0;
                if (a == 1) return r;

                if ((a % 4 == 3) && (n % 4) == 3) r = -r;
                std::swap(a, n);
            }
        }

        // Performs a strong lucas probable primality test on n. Assumes n is odd and has no tiny
        // factors.
        inline bool strong_lucas_probable_prime(uint64_t n, uint64_t mont1, uint64_t nneginv) {
            int64_t Dc = 5; // Candidate for D.
            int j = jacobi(Dc, n);

            while (j != -1) {
                if (j == 0) return false;
                if (Dc < 0) Dc =  2 - Dc;
                else        Dc = -2 - Dc;

                // Check if n is a square if we don't find a proper D quickly.
                if (Dc == 13) {
                    uint64_t n_sqrt = op::isqrt(n);
                    if (n_sqrt * n_sqrt == n) return false;
                }

                j = jacobi(Dc, n);
            }

            // We need to divide by two, convert modular inverse of 2 to Montgomery form.
            const uint64_t mont_div_2 = op::divu128_64((n + 1) / 2, 0, n).second;
            
            // Set up primality test.
            uint64_t D = op::divu128_64(Dc < 0 ? n + Dc : Dc, 0, n).second;
            uint64_t Q = op::divu128_64(Dc > 1 ? n + (1 - Dc) / 4 : (1 - Dc) / 4, 0, n).second;
            uint64_t U = mont1;
            uint64_t V = mont1;
            uint64_t QQ = Q;
            
            uint64_t d = n + 1; // Can't overflow because 2^64-1 is divisible by 3.
            int s = 0;
            while ((d & 1) == 0) {
                s += 1;
                d >>= 1;
            }

            // Mask to step through bits from left to right.
            uint64_t mask = d;
            mask |= mask >> 1; mask |= mask >> 2; mask |= mask >> 4;
            mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;
            mask = (mask >> 1) ^ (mask >> 2);

            while (mask) {
                U = detail::montmul64(U, V, n, nneginv);
                V = detail::montmul64(V, V, n, nneginv);
                V = detail::submod64(V, QQ, n);
                V = detail::submod64(V, QQ, n);
                QQ = detail::montmul64(QQ, QQ, n, nneginv);

                if (d & mask) {
                    uint64_t U2 = detail::addmod64(U, V, n);
                    uint64_t V2 = detail::addmod64(detail::montmul64(D, U, n, nneginv), V, n);
                    U = detail::montmul64(U2, mont_div_2, n, nneginv);
                    V = detail::montmul64(V2, mont_div_2, n, nneginv);
                    QQ = detail::montmul64(QQ, Q, n, nneginv);
                }

                mask >>= 1;
            }
            
            if (U == 0 || V == 0) return true;

            for (int r = 1; r < s; ++r) {
                V = detail::montmul64(V, V, n, nneginv);
                V = detail::submod64(V, QQ, n);
                V = detail::submod64(V, QQ, n);
                
                if (V == 0) return true;

                QQ = detail::montmul64(QQ, QQ, n, nneginv);
            }

            return false;
        }
    }


    template<class T>
    inline bool is_prime(T n) {
        static_assert(std::is_integral<T>::value, "op::is_prime only works on integer types");

        if (n < 2) return false;
        return op::is_prime<uint64_t>(n);
    }


    template<>
    inline bool is_prime(uint64_t n) {
        constexpr const int max_smallprimeset = 10000; // About a 4K table.
        static std::set<unsigned> smallprimeset;

        if (smallprimeset.size() == 0) {
            op::primes_below(max_smallprimeset,
                             std::inserter(smallprimeset, smallprimeset.begin()));
        }

        if (n <= 3) return n >= 2;
        if ((n & 1) == 0) return false;
        if (n < max_smallprimeset) return smallprimeset.count(n);

        // Check divisibility by small odd primes.
        if (n * 0xaaaaaaaaaaaaaaabull <= 0x5555555555555555ull) return false; //  3 
        if (n * 0xcccccccccccccccdull <= 0x3333333333333333ull) return false; //  5 
        if (n * 0x6db6db6db6db6db7ull <= 0x2492492492492492ull) return false; //  7 
        if (n * 0x2e8ba2e8ba2e8ba3ull <= 0x1745d1745d1745d1ull) return false; // 11 
        if (n * 0x4ec4ec4ec4ec4ec5ull <= 0x13b13b13b13b13b1ull) return false; // 13 
        if (n * 0xf0f0f0f0f0f0f0f1ull <= 0x0f0f0f0f0f0f0f0full) return false; // 17 
        if (n * 0x86bca1af286bca1bull <= 0x0d79435e50d79435ull) return false; // 19 
        if (n * 0xd37a6f4de9bd37a7ull <= 0x0b21642c8590b216ull) return false; // 23 
        if (n * 0x34f72c234f72c235ull <= 0x08d3dcb08d3dcb08ull) return false; // 29 
        if (n * 0xef7bdef7bdef7bdfull <= 0x0842108421084210ull) return false; // 31 
        if (n * 0x14c1bacf914c1badull <= 0x06eb3e45306eb3e4ull) return false; // 37 
        if (n * 0x8f9c18f9c18f9c19ull <= 0x063e7063e7063e70ull) return false; // 41 
        if (n * 0x82fa0be82fa0be83ull <= 0x05f417d05f417d05ull) return false; // 43 
        if (n * 0x51b3bea3677d46cfull <= 0x0572620ae4c415c9ull) return false; // 47 

        // Bring out the big guns (deterministic Miller-Rabin or strong lucas probable prime check).
        uint64_t mont1 = op::divu128_64(1, 0, n).second;
        uint64_t montn1 = op::divu128_64(n-1, 0, n).second;
        uint64_t nneginv = detail::mont_modinv(n).second;

        if (n < 341531ull) {
            return detail::strong_probable_prime(n, 9345883071009581737ull, mont1, montn1, nneginv);
        }

        if (n < 1050535501ull) {
            return
                detail::strong_probable_prime(n,     336781006125ull, mont1, montn1, nneginv) &&
                detail::strong_probable_prime(n, 9639812373923155ull, mont1, montn1, nneginv);
        }

        if (n < 350269456337ull) {
            return
                detail::strong_probable_prime(n,  4230279247111683200ull, mont1, montn1, nneginv) &&
                detail::strong_probable_prime(n, 14694767155120705706ull, mont1, montn1, nneginv) &&
                detail::strong_probable_prime(n, 16641139526367750375ull, mont1, montn1, nneginv);
        }

        // Once we would need to check 4 bases a strong lucas probable prime check is faster.
        return detail::strong_probable_prime(n, 2, mont1, montn1, nneginv) &&
               detail::strong_lucas_probable_prime(n, mont1, nneginv);
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
            uint64_t nneginv = detail::mont_modinv(n).second;
            uint64_t mont1 = op::divu128_64(1, 0, n).second;

            uint64_t g, r, q, x, ys;
            uint64_t rng = 42;

            do {
                uint64_t y = rng++;
                uint64_t m = 256;

                g = q = r = 1;
                q = mont1;

                do {
                    x = y;
                    for (uint64_t i = 0; i < r; ++i) {
                        // Checking for a potential overflow here is slower than to just waste a
                        // trial when it overflows.
                        y = detail::montmul64(y, y, n, nneginv) + 1;
                    }

                    for (uint64_t k = 0; k < r && g == 1; k += m) {
                        ys = y;
                        for (uint64_t i = 0; q && i < std::min(m, r-k); ++i) {
                            y = detail::montmul64(y, y, n, nneginv) + 1;
                            q = detail::montmul64(q, x < y ? y-x : x-y, n, nneginv);
                        }

                        g = op::gcd(q, n);
                    }

                    r *= 2;
                } while (g == 1);

                if (g == n) {
                    do {
                        ys = detail::montmul64(ys, ys, n, nneginv) + 1;
                        g = op::gcd(x < ys ? ys-x : x-ys, n);
                    } while (g == 1);
                }
            } while (g == n);
  
            return g;
        }
    }


    inline std::vector<uint64_t> prime_factors(uint64_t n) {
        static const std::array<std::tuple<int, uint64_t, uint64_t>, 167> smallprimes {
            std::make_tuple(  3, 0xaaaaaaaaaaaaaaabull, 0x5555555555555555ull),
            std::make_tuple(  5, 0xcccccccccccccccdull, 0x3333333333333333ull),
            std::make_tuple(  7, 0x6db6db6db6db6db7ull, 0x2492492492492492ull),
            std::make_tuple( 11, 0x2e8ba2e8ba2e8ba3ull, 0x1745d1745d1745d1ull),
            std::make_tuple( 13, 0x4ec4ec4ec4ec4ec5ull, 0x13b13b13b13b13b1ull),
            std::make_tuple( 17, 0xf0f0f0f0f0f0f0f1ull, 0x0f0f0f0f0f0f0f0full),
            std::make_tuple( 19, 0x86bca1af286bca1bull, 0x0d79435e50d79435ull),
            std::make_tuple( 23, 0xd37a6f4de9bd37a7ull, 0x0b21642c8590b216ull),
            std::make_tuple( 29, 0x34f72c234f72c235ull, 0x08d3dcb08d3dcb08ull),
            std::make_tuple( 31, 0xef7bdef7bdef7bdfull, 0x0842108421084210ull),
            std::make_tuple( 37, 0x14c1bacf914c1badull, 0x06eb3e45306eb3e4ull),
            std::make_tuple( 41, 0x8f9c18f9c18f9c19ull, 0x063e7063e7063e70ull),
            std::make_tuple( 43, 0x82fa0be82fa0be83ull, 0x05f417d05f417d05ull),
            std::make_tuple( 47, 0x51b3bea3677d46cfull, 0x0572620ae4c415c9ull),
            std::make_tuple( 53, 0x21cfb2b78c13521dull, 0x04d4873ecade304dull),
            std::make_tuple( 59, 0xcbeea4e1a08ad8f3ull, 0x0456c797dd49c341ull),
            std::make_tuple( 61, 0x4fbcda3ac10c9715ull, 0x04325c53ef368eb0ull),
            std::make_tuple( 67, 0xf0b7672a07a44c6bull, 0x03d226357e16ece5ull),
            std::make_tuple( 71, 0x193d4bb7e327a977ull, 0x039b0ad12073615aull),
            std::make_tuple( 73, 0x7e3f1f8fc7e3f1f9ull, 0x0381c0e070381c0eull),
            std::make_tuple( 79, 0x9b8b577e613716afull, 0x033d91d2a2067b23ull),
            std::make_tuple( 83, 0xa3784a062b2e43dbull, 0x03159721ed7e7534ull),
            std::make_tuple( 89, 0xf47e8fd1fa3f47e9ull, 0x02e05c0b81702e05ull),
            std::make_tuple( 97, 0xa3a0fd5c5f02a3a1ull, 0x02a3a0fd5c5f02a3ull),
            std::make_tuple(101, 0x3a4c0a237c32b16dull, 0x0288df0cac5b3f5dull),
            std::make_tuple(103, 0xdab7ec1dd3431b57ull, 0x027c45979c95204full),
            std::make_tuple(107, 0x77a04c8f8d28ac43ull, 0x02647c69456217ecull),
            std::make_tuple(109, 0xa6c0964fda6c0965ull, 0x02593f69b02593f6ull),
            std::make_tuple(113, 0x90fdbc090fdbc091ull, 0x0243f6f0243f6f02ull),
            std::make_tuple(127, 0x7efdfbf7efdfbf7full, 0x0204081020408102ull),
            std::make_tuple(131, 0x03e88cb3c9484e2bull, 0x01f44659e4a42715ull),
            std::make_tuple(137, 0xe21a291c077975b9ull, 0x01de5d6e3f8868a4ull),
            std::make_tuple(139, 0x3aef6ca970586723ull, 0x01d77b654b82c339ull),
            std::make_tuple(149, 0xdf5b0f768ce2cabdull, 0x01b7d6c3dda338b2ull),
            std::make_tuple(151, 0x6fe4dfc9bf937f27ull, 0x01b2036406c80d90ull),
            std::make_tuple(157, 0x5b4fe5e92c0685b5ull, 0x01a16d3f97a4b01aull),
            std::make_tuple(163, 0x1f693a1c451ab30bull, 0x01920fb49d0e228dull),
            std::make_tuple(167, 0x8d07aa27db35a717ull, 0x01886e5f0abb0499ull),
            std::make_tuple(173, 0x882383b30d516325ull, 0x017ad2208e0ecc35ull),
            std::make_tuple(179, 0xed6866f8d962ae7bull, 0x016e1f76b4337c6cull),
            std::make_tuple(181, 0x3454dca410f8ed9dull, 0x016a13cd15372904ull),
            std::make_tuple(191, 0x1d7ca632ee936f3full, 0x01571ed3c506b39aull),
            std::make_tuple(193, 0x70bf015390948f41ull, 0x015390948f40feacull),
            std::make_tuple(197, 0xc96bdb9d3d137e0dull, 0x014cab88725af6e7ull),
            std::make_tuple(199, 0x2697cc8aef46c0f7ull, 0x0149539e3b2d066eull),
            std::make_tuple(211, 0xc0e8f2a76e68575bull, 0x013698df3de07479ull),
            std::make_tuple(223, 0x687763dfdb43bb1full, 0x0125e22708092f11ull),
            std::make_tuple(227, 0x1b10ea929ba144cbull, 0x0120b470c67c0d88ull),
            std::make_tuple(229, 0x1d10c4c0478bbcedull, 0x011e2ef3b3fb8744ull),
            std::make_tuple(233, 0x63fb9aeb1fdcd759ull, 0x0119453808ca29c0ull),
            std::make_tuple(239, 0x64afaa4f437b2e0full, 0x0112358e75d30336ull),
            std::make_tuple(241, 0xf010fef010fef011ull, 0x010fef010fef010full),
            std::make_tuple(251, 0x28cbfbeb9a020a33ull, 0x0105197f7d734041ull),
            std::make_tuple(257, 0xff00ff00ff00ff01ull, 0x00ff00ff00ff00ffull),
            std::make_tuple(263, 0xd624fd1470e99cb7ull, 0x00f92fb2211855a8ull),
            std::make_tuple(269, 0x8fb3ddbd6205b5c5ull, 0x00f3a0d52cba8723ull),
            std::make_tuple(271, 0xd57da36ca27acdefull, 0x00f1d48bcee0d399ull),
            std::make_tuple(277, 0xee70c03b25e4463dull, 0x00ec979118f3fc4dull),
            std::make_tuple(281, 0xc5b1a6b80749cb29ull, 0x00e939651fe2d8d3ull),
            std::make_tuple(283, 0x47768073c9b97113ull, 0x00e79372e225fe30ull),
            std::make_tuple(293, 0x2591e94884ce32adull, 0x00dfac1f74346c57ull),
            std::make_tuple(307, 0xf02806abc74be1fbull, 0x00d578e97c3f5fe5ull),
            std::make_tuple(311, 0x7ec3e8f3a7198487ull, 0x00d2ba083b445250ull),
            std::make_tuple(313, 0x58550f8a39409d09ull, 0x00d161543e28e502ull),
            std::make_tuple(317, 0xec9e48ae6f71de15ull, 0x00cebcf8bb5b4169ull),
            std::make_tuple(331, 0x2ff3a018bfce8063ull, 0x00c5fe740317f9d0ull),
            std::make_tuple(337, 0x7f9ec3fcf61fe7b1ull, 0x00c2780613c0309eull),
            std::make_tuple(347, 0x89f5abe570e046d3ull, 0x00bcdd535db1cc5bull),
            std::make_tuple(349, 0xda971b23f1545af5ull, 0x00bbc8408cd63069ull),
            std::make_tuple(353, 0x79d5f00b9a7862a1ull, 0x00b9a7862a0ff465ull),
            std::make_tuple(359, 0x4dba1df32a128a57ull, 0x00b68d31340e4307ull),
            std::make_tuple(367, 0x87530217b7747d8full, 0x00b2927c29da5519ull),
            std::make_tuple(373, 0x30baae53bb5e06ddull, 0x00afb321a1496fdfull),
            std::make_tuple(379, 0xee70206c12e9b5b3ull, 0x00aceb0f891e6551ull),
            std::make_tuple(383, 0xcdde9462ec9dbe7full, 0x00ab1cbdd3e2970full),
            std::make_tuple(389, 0xafb64b05ec41cf4dull, 0x00a87917088e262bull),
            std::make_tuple(397, 0x02944ff5aec02945ull, 0x00a513fd6bb00a51ull),
            std::make_tuple(401, 0x2cb033128382df71ull, 0x00a36e71a2cb0331ull),
            std::make_tuple(409, 0x1ccacc0c84b1c2a9ull, 0x00a03c1688732b30ull),
            std::make_tuple(419, 0x19a93db575eb3a0bull, 0x009c69169b30446dull),
            std::make_tuple(421, 0xcebeef94fa86fe2dull, 0x009baade8e4a2f6eull),
            std::make_tuple(431, 0x6faa77fb3f8df54full, 0x00980e4156201301ull),
            std::make_tuple(433, 0x68a58af00975a751ull, 0x00975a750ff68a58ull),
            std::make_tuple(439, 0xd56e36d0c3efac07ull, 0x009548e4979e0829ull),
            std::make_tuple(443, 0xd8b44c47a8299b73ull, 0x0093efd1c50e726bull),
            std::make_tuple(449, 0x02d9ccaf9ba70e41ull, 0x0091f5bcb8bb02d9ull),
            std::make_tuple(457, 0x0985e1c023d9e879ull, 0x008f67a1e3fdc261ull),
            std::make_tuple(461, 0x2a343316c494d305ull, 0x008e2917e0e702c6ull),
            std::make_tuple(463, 0x70cb7916ab67652full, 0x008d8be33f95d715ull),
            std::make_tuple(467, 0xd398f132fb10fe5bull, 0x008c55841c815ed5ull),
            std::make_tuple(479, 0x6f2a38a6bf54fa1full, 0x0088d180cd3a4133ull),
            std::make_tuple(487, 0x211df689b98f81d7ull, 0x00869222b1acf1ceull),
            std::make_tuple(491, 0x0e994983e90f1ec3ull, 0x0085797b917765abull),
            std::make_tuple(499, 0xad671e44bed87f3bull, 0x008355ace3c897dbull),
            std::make_tuple(503, 0xf9623a0516e70fc7ull, 0x00824a4e60b3262bull),
            std::make_tuple(509, 0x4b7129be9dece355ull, 0x0080c121b28bd1baull),
            std::make_tuple(521, 0x190f3b7473f62c39ull, 0x007dc9f3397d4c29ull),
            std::make_tuple(523, 0x63dacc9aad46f9a3ull, 0x007d4ece8fe88139ull),
            std::make_tuple(541, 0xc1108fda24e8d035ull, 0x0079237d65bcce50ull),
            std::make_tuple(547, 0xb77578472319bd8bull, 0x0077cf53c5f7936cull),
            std::make_tuple(557, 0x473d20a1c7ed9da5ull, 0x0075a8accfbdd11eull),
            std::make_tuple(563, 0xfbe85af0fea2c8fbull, 0x007467ac557c228eull),
            std::make_tuple(569, 0x58a1f7e6ce0f4c09ull, 0x00732d70ed8db8e9ull),
            std::make_tuple(571, 0x1a00e58c544986f3ull, 0x0072c62a24c3797full),
            std::make_tuple(577, 0x7194a17f55a10dc1ull, 0x007194a17f55a10dull),
            std::make_tuple(587, 0x7084944785e33763ull, 0x006fa549b41da7e7ull),
            std::make_tuple(593, 0xba10679bd84886b1ull, 0x006e8419e6f61221ull),
            std::make_tuple(599, 0xebe9c6bb31260967ull, 0x006d68b5356c207bull),
            std::make_tuple(601, 0x97a3fe4bd1ff25e9ull, 0x006d0b803685c01bull),
            std::make_tuple(607, 0x6c6388395b84d99full, 0x006bf790a8b2d207ull),
            std::make_tuple(613, 0x8c51da6a1335df6dull, 0x006ae907ef4b96c2ull),
            std::make_tuple(617, 0x46f3234475d5add9ull, 0x006a37991a23aeadull),
            std::make_tuple(619, 0x905605ca3c619a43ull, 0x0069dfbdd4295b66ull),
            std::make_tuple(631, 0xcee8dff304767747ull, 0x0067dc4c45c8033eull),
            std::make_tuple(641, 0xff99c27f00663d81ull, 0x00663d80ff99c27full),
            std::make_tuple(643, 0xacca407f671ddc2bull, 0x0065ec17e3559948ull),
            std::make_tuple(647, 0xe71298bac1e12337ull, 0x00654ac835cfba5cull),
            std::make_tuple(653, 0xfa1e94309cd09045ull, 0x00645c854ae10772ull),
            std::make_tuple(659, 0xbebccb8e91496b9bull, 0x006372990e5f901full),
            std::make_tuple(661, 0x312fa30cc7d7b8bdull, 0x006325913c07beefull),
            std::make_tuple(673, 0x6160ff9e9f006161ull, 0x006160ff9e9f0061ull),
            std::make_tuple(677, 0x6b03673b5e28152dull, 0x0060cdb520e5e88eull),
            std::make_tuple(683, 0xfe802ffa00bfe803ull, 0x005ff4017fd005ffull),
            std::make_tuple(691, 0xe66fe25c9e907c7bull, 0x005ed79e31a4dccdull),
            std::make_tuple(701, 0x3f8b236c76528895ull, 0x005d7d42d48ac5efull),
            std::make_tuple(709, 0xf6f923bf01ce2c0dull, 0x005c6f35ccba5028ull),
            std::make_tuple(719, 0x6c3d3d98bed7c42full, 0x005b2618ec6ad0a5ull),
            std::make_tuple(727, 0x30981efcd4b010e7ull, 0x005a2553748e42e7ull),
            std::make_tuple(733, 0x6f691fc81ebbe575ull, 0x0059686cf744cd5bull),
            std::make_tuple(739, 0xb10480ddb47b52cbull, 0x0058ae97bab79976ull),
            std::make_tuple(743, 0x74cd59ed64f3f0d7ull, 0x0058345f1876865full),
            std::make_tuple(751, 0x0105cb81316d6c0full, 0x005743d5bb24795aull),
            std::make_tuple(757, 0x9be64c6d91c1195dull, 0x005692c4d1ab74abull),
            std::make_tuple(761, 0x71b3f945a27b1f49ull, 0x00561e46a4d5f337ull),
            std::make_tuple(769, 0x77d80d50e508fd01ull, 0x005538ed06533997ull),
            std::make_tuple(773, 0xa5eb778e133551cdull, 0x0054c807f2c0bec2ull),
            std::make_tuple(787, 0x18657d3c2d8a3f1bull, 0x005345efbc572d36ull),
            std::make_tuple(797, 0x2e40e220c34ad735ull, 0x00523a758f941345ull),
            std::make_tuple(809, 0xa76593c70a714919ull, 0x005102370f816c89ull),
            std::make_tuple(811, 0x1eef452124eea383ull, 0x0050cf129fb94acfull),
            std::make_tuple(821, 0x38206dc242ba771dull, 0x004fd31941cafdd1ull),
            std::make_tuple(823, 0x4cd4c35807772287ull, 0x004fa1704aa75945ull),
            std::make_tuple(827, 0x83de917d5e69ddf3ull, 0x004f3ed6d45a63adull),
            std::make_tuple(829, 0x882ef0403b4a6c15ull, 0x004f0de57154ebedull),
            std::make_tuple(839, 0xf8fb6c51c606b677ull, 0x004e1cae8815f811ull),
            std::make_tuple(853, 0xb4abaac446d3e1fdull, 0x004cd47ba5f6ff19ull),
            std::make_tuple(857, 0xa9f83bbe484a14e9ull, 0x004c78ae734df709ull),
            std::make_tuple(859, 0x0bebbc0d1ce874d3ull, 0x004c4b19ed85cfb8ull),
            std::make_tuple(863, 0xbd418eaf0473189full, 0x004bf093221d1218ull),
            std::make_tuple(877, 0x44e3af6f372b7e65ull, 0x004aba3c21dc633full),
            std::make_tuple(881, 0xc87fdace4f9e5d91ull, 0x004a6360c344de00ull),
            std::make_tuple(883, 0xec93479c446bd9bbull, 0x004a383e9f74d68aull),
            std::make_tuple(887, 0xdac4d592e777c647ull, 0x0049e28fbabb9940ull),
            std::make_tuple(907, 0xa63ea8c8f61f0c23ull, 0x0048417b57c78cd7ull),
            std::make_tuple(911, 0xe476062ea5cbbb6full, 0x0047f043713f3a2bull),
            std::make_tuple(919, 0xdf68761c69daac27ull, 0x00474ff2a10281cfull),
            std::make_tuple(929, 0xb813d737637aa061ull, 0x00468b6f9a978f91ull),
            std::make_tuple(937, 0xa3a77aac1fb15099ull, 0x0045f13f1caff2e2ull),
            std::make_tuple(941, 0x17f0c3e0712c5825ull, 0x0045a5228cec23e9ull),
            std::make_tuple(947, 0xfd912a70ff30637bull, 0x0045342c556c66b9ull),
            std::make_tuple(953, 0xfbb3b5dc01131289ull, 0x0044c4a23feeced7ull),
            std::make_tuple(967, 0x856d560a0f5acdf7ull, 0x0043c5c20d3c9fe6ull),
            std::make_tuple(971, 0x96472f314d3f89e3ull, 0x00437e494b239798ull),
            std::make_tuple(977, 0xa76f5c7ed2253531ull, 0x0043142d118e47cbull),
            std::make_tuple(983, 0x816eae7c7bf69fe7ull, 0x0042ab5c73a13458ull),
            std::make_tuple(991, 0xb6a2bea4cfb1781full, 0x004221950db0f3dbull),
            std::make_tuple(997, 0xa3900c53318e81edull, 0x0041bbb2f80a4553ull),
        };

        if (n <= 1) return {1};

        std::vector<uint64_t> factors;
        while ((n & 1) == 0) {
            factors.push_back(2);
            n >>= 1;
        }

        for (auto checker : smallprimes) {
            while (n * std::get<1>(checker) <= std::get<2>(checker)) {
                factors.push_back(std::get<0>(checker));
                n *= std::get<1>(checker);
            }

            if (uint64_t(std::get<0>(checker)) > n) break;
        }

        std::vector<uint64_t> to_factor = {n};
        while (to_factor.size()) {
            n = to_factor.back();
            to_factor.pop_back();

            if (n == 1) continue;
            if (op::is_prime(n)) {
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


    inline std::map<uint64_t, int> factorization(uint64_t n) {
        std::map<uint64_t, int> factors;

        for (auto factor : op::prime_factors(n)) factors[factor]++;

        return factors;
    }
    

    namespace detail {
        template<class T, class U,
                 bool=(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value),
                 bool=std::is_signed<T>::value, bool=std::is_signed<U>::value,
                 bool=std::is_floating_point<T>::value, bool=std::is_floating_point<U>::value>
        struct safe_less_helper;

        // One of the types not arithmetic - just use regular comparison.
        template<class T, class U, bool a, bool b, bool c, bool d>
        struct safe_less_helper<T, U, false, a, b, c, d> {
            static constexpr bool eval(const T& x, const U& y) { return x < y; }
        };

        template<class T, class U, bool same_sign, bool same_is_float>
        struct safe_less_helper<T, U, true, same_sign, same_sign, same_is_float, same_is_float> {
            static constexpr bool eval(const T& x, const U& y) { return x < y; }
        };

        template<class T, class U>
        struct safe_less_helper<T, U, true, true, false, false, false> {
            static constexpr bool eval(const T& x, const U& y) {
                return x < 0 || static_cast<typename std::make_unsigned<T>::type>(x) < y;
            }
        };

        template<class T, class U>
        struct safe_less_helper<T, U, true, false, true, false, false> {
            static constexpr bool eval(const T& x, const U& y) {
                return y > 0 && x < static_cast<typename std::make_unsigned<U>::type>(y);
            }
        };

        template<class T, class U, bool a, bool b>
        struct safe_less_helper<T, U, true, a, b, true, false> {
            static constexpr bool eval(const T& x, const U& y) {
                return
                    // This should never give a false positive. A false positive occurs when the
                    // real value of y >= x, but this comparison returns true. If y > x then an
                    // adjacent floating point number representation could be x, but not a number
                    // smaller than x.  If y == x the conversion would be exact.
                    x < T(y) ?                                            
                        true : (
                    // Similar story here, a false positive is impossible.
                    x > T(y) ?
                        false : (
                    // We can safely convert to int.
                    T(std::numeric_limits<U>::min()) < x && x < T(std::numeric_limits<U>::max()) ?
                        U(x) < y : (
                    // NaN comparison is always false.
                    std::isnan(x) ?
                        false : 
                    // Now it's tricky. The floating point value could be slightly under or slightly
                    // above std::numerical_limits<U>::max or min. After division by the radix it
                    // should be safe to exactly compare as integers.
                        is_less_div_radix<std::numeric_limits<T>::radix>(x, y))));
            }

            template<int r>
            static constexpr bool is_less_div_radix(const T& x, const U& y) {
                return U(x/r) < y/r || U(x/r) == y/r && U((x/r-U(x/r))*r) < y%r;
            }
        };

        template<class T, class U, bool a, bool b>
        struct safe_less_helper<T, U, true, a, b, false, true> {
            static constexpr bool eval(const T& x, const U& y) {
                return
                    U(x) < y ?                                            
                        true : (
                    U(x) > y ?
                        false : (
                    U(std::numeric_limits<T>::min()) < y && y < U(std::numeric_limits<T>::max()) ?
                        x < T(y) : (
                    std::isnan(y) ?
                        false : 
                        is_less_div_radix<std::numeric_limits<U>::radix>(x, y))));
            }

            template<int r>
            static constexpr bool is_less_div_radix(const T& x, const U& y) {
                return x/r < T(y/r) || x/r == T(y/r) && x%r < T((y/r-T(y/r))*r);
            }
        };
    }

    template<class T, class U>
    inline constexpr bool safe_less(const T& x, const U& y) {
        return detail::safe_less_helper<T, U>::eval(x, y);
    }

    template<class T, class U>
    inline constexpr bool safe_greater(const T& x, const U& y) { return safe_less(y, x); }

    template<class T, class U>
    inline constexpr bool safe_less_equal(const T& x, const U& y) { return !safe_less(y, x); }

    template<class T, class U>
    inline constexpr bool safe_greater_equal(const T& x, const U& y) { return !safe_less(x, y); }
}

#endif
