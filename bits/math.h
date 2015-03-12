#ifndef OP_MATH_H
#define OP_MATH_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <set>
#include <iterator>
#include <map>

#include "intrin.h"


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
    bool isprime(T n);

    // Returns a vector containing the prime factors of n.
    std::vector<uint64_t> prime_factors(uint64_t n);

    // Returns a map with the prime factors of n as the key, and how many times the prime factor
    // occurs in the factorization of n as the value.
    std::map<uint64_t, int> factorization(uint64_t n);    

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


    int64_t ipow(int32_t base, uint8_t exp) {
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


        // Computes aR * bR mod N with R = 2**64.
        inline uint64_t montmul64(uint64_t a, uint64_t b, uint64_t N, uint64_t Nneginv) {
            uint64_t Th, Tl, m, mNh, mNl, th;

            std::tie(Th, Tl) = op::mulu64(a, b);
            m = Tl * Nneginv;
            std::tie(mNh, mNl) = op::mulu64(m, N);

            bool lc = Tl + mNl < Tl;
            th = Th + mNh + lc;
            bool hc = (th < Th) || (th == Th && lc);

            if (hc > 0 || th >= N) th = th - N;

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
            uint64_t b = op::modu128_64(base % n, 0, n);
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
            const uint64_t mont_div_2 = op::modu128_64((n + 1) / 2, 0, n);
            
            // Set up primality test.
            uint64_t D = op::modu128_64(Dc < 0 ? n + Dc : Dc, 0, n);
            uint64_t Q = op::modu128_64(Dc > 1 ? n + (1 - Dc) / 4 : (1 - Dc) / 4, 0, n);
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
    inline bool isprime(T n) {
        static_assert(std::is_integral<T>::value, "op::isprime only works on integer types");

        if (n < 2) return false;
        return op::isprime<uint64_t>(n);
    }


    template<>
    inline bool isprime(uint64_t n) {
        constexpr int max_smallprimeset = 10000; // About a 4K table.
        static std::set<unsigned> smallprimeset;

        if (smallprimeset.size() == 0) {
            op::primes_below(max_smallprimeset, std::inserter(smallprimeset, smallprimeset.begin()));
        }

        if (n <= 3) return n >= 2;
        if (n % 2 == 0) return false;
        if (n < max_smallprimeset) return smallprimeset.count(n);

        // Looping this is significantly slower.
        if (n % 3 == 0 || n % 5 == 0 || n % 7 == 0 || n % 11 == 0 || n % 13 == 0 || n % 17 == 0 ||
            n % 19 == 0 || n % 23 == 0 || n % 29 == 0 || n % 31 == 0 || n % 47 == 0) return false;

        uint64_t mont1 = op::modu128_64(1, 0, n);
        uint64_t montn1 = op::modu128_64(n-1, 0, n);
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
            uint64_t mont1 = op::modu128_64(1, 0, n);

            uint64_t g, r, q, x, ys;
            uint64_t rng = 1;

            do {
                uint64_t y = rng++;
                uint64_t m = 100;

                g = q = r = 1;
                q = mont1;

                do {
                    x = y;
                    for (uint64_t i = 0; i < r; ++i) {
                        y = detail::montmul64(y, y, n, nneginv) + 1;
                    }

                    for (uint64_t k = 0; k < r && g == 1; k += m) {
                        ys = y;
                        for (uint64_t i = 0; i < std::min(m, r-k); ++i) {
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
        static std::vector<int> smallprimes;
        if (smallprimes.size() == 0) op::primes_below(1000, std::back_inserter(smallprimes));

        if (n <= 1) return {1};

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


    inline std::map<uint64_t, int> factorization(uint64_t n) {
        std::map<uint64_t, int> factors;

        for (auto factor : op::prime_factors(n)) factors[factor]++;

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
