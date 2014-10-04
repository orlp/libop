#ifndef OP_UTILITY_H
#define OP_UTILITY_H

#include <array>
#include <type_traits>
#include <utility>


namespace op {
    namespace detail {
        template<std::size_t ...> struct indices {};

        template<std::size_t N, std::size_t ...I>
        struct build_indices : build_indices<N-1, N-1, I...> { };

        template<std::size_t ...I>
        struct build_indices<0, I...> : indices<I...> { };

        template<class T, std::size_t N, std::size_t... i>
        inline constexpr std::array<T, N> array_unpacker(T (&a)[N], indices<i...>) {
            return {{ a[i]... }};
        }
    }

    // make_array(...) is a helper function to create a std::array with a derived type and size.
    template<class... Types>
    inline constexpr std::array<typename std::common_type<Types...>::type, sizeof...(Types)>
    make_array(Types&&... t) {
        return {{ std::forward<Types>(t)... }};
    }


    // make_array<CT>(...) is a helper function to create a std::array with a derived size.
    template<class CT, class... Types>
    inline constexpr std::array<CT, sizeof...(Types)> make_array(Types&&... t) {
        return {{ static_cast<CT>(std::forward<Types>(t))... }};
    }


    // to_array(c_array) is a helper function to create a std::array from a C array.
    template<class T, std::size_t N>
    inline constexpr std::array<T, N> to_array(T (&a)[N]) {
        return detail::array_unpacker(a, detail::build_indices<N>());
    }
}


#endif
