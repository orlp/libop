#ifndef OP_UTILITY_H
#define OP_UTILITY_H

#include <type_traits>
#include <utility>
#include <array>

namespace op {
    template<class... T>
    constexpr std::array<typename std::common_type<T...>::type, sizeof...(T)> make_array(T&&... t) {
        return {{ std::forward<T>(t)... }};
    }

    template<class CT, class... T>
    constexpr std::array<CT, sizeof...(T)> make_array(T&&... t) {
        return {{ static_cast<CT>(std::forward<T>(t))... }};
    }
}


#endif