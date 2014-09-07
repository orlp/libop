#ifndef OP_ILOG_H
#define OP_ILOG_H

#include <limits>

namespace op {
    namespace detail {
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

        template<class T, T base>
        inline constexpr int max_pot_power(int result = 1) {
            return (base <= (std::numeric_limits<T>::max() / base)) ? max_pot_power<T, base * base>(result * 2) : result; 
        }
    }

    template<int exp, class T>
    inline constexpr T pow(T base) {
        static_assert(exp < 64, "pow exponents >= 64 can only overflow");

        return exp < 0 ? 1 / detail::pow_impl<-exp>(base) : detail::pow_impl<exp>(base);
    }




    template<int base, class T>
    inline  T ilog(T n) {
        // int i = -(n == 0);

        int x = 0;

        for (int i = detail::max_pot_power<T, base>(); i; i /= 2) x ++;

        // for (constexpr T cmp = detail::max_power<T>(base); cmp; cmp /= base) {
        //     if (n >= cmp) {
        //         i += k;
        //         n /= cmp;
        //     }
        // }

        return x;
    }
}

#endif