#ifndef OP_UTILITY_H
#define OP_UTILITY_H

#include <array>
#include <type_traits>
#include <utility>

namespace op {
    template<class... T>
    inline constexpr std::array<typename std::common_type<T...>::type, sizeof...(T)> make_array(T&&... t) {
        return {{ std::forward<T>(t)... }};
    }

    template<class CT, class... T>
    inline constexpr std::array<CT, sizeof...(T)> make_array(T&&... t) {
        return {{ static_cast<CT>(std::forward<T>(t))... }};
    }
}


#endif
