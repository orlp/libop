#ifndef OP_TYPE_TRAITS_H
#define OP_TYPE_TRAITS_H

#include <cstddef>
#include <type_traits>


namespace op {
    // Inherits from std::true_type if T and U are Swappable, std::false_type if otherwise.
    template<class T, class U = T>
    struct is_swappable;

    // Same as is_swappable, with the added condition that the swap must be noexcept.
    template<class T, class U = T>
    struct is_nothrow_swappable;
}



// Implementation.
namespace op {
    namespace detail {
        namespace swap_adl_tests {
            // If swap ADL finds this then it would call std::swap otherwise (same signature).
            struct tag {};

            template<class T> tag swap(T&, T&) noexcept(
                std::is_nothrow_move_constructible<T>::value &&
                std::is_nothrow_move_assignable<T>::value
            );

            template<class T, std::size_t N>
            tag swap(T (&a)[N], T (&b)[N]) noexcept(noexcept(swap(*a, *b)));

            // Helper functions to test if an unqualified swap is possible, and whether it becomes a
            // std::swap.
            template<class, class> std::false_type uses_std(...);
            template<class T, class U>
            std::is_same<decltype(swap(std::declval<T&>(), std::declval<U&>())), tag> uses_std(int);

            template<class, class> std::false_type can_swap(...) noexcept(false);
            template<class T, class U,
                     class = decltype(swap(std::declval<T&>(), std::declval<U&>()))>
            std::true_type can_swap(int)
                    noexcept(noexcept(swap(std::declval<T&>(), std::declval<U&>())));
        }
    }

    template<class T, class U>
    struct is_swappable : std::integral_constant<bool, 
        decltype(detail::swap_adl_tests::can_swap<T, U>(0))::value &&
            (!decltype(detail::swap_adl_tests::uses_std<T, U>(0))::value ||
                (std::is_move_assignable<T>::value && std::is_move_constructible<T>::value))
    > {};

    template<class T, std::size_t N>
    struct is_swappable<T[N], T[N]> : std::integral_constant<bool, 
        decltype(detail::swap_adl_tests::can_swap<T[N], T[N]>(0))::value &&
            (!decltype(detail::swap_adl_tests::uses_std<T[N], T[N]>(0))::value ||
                is_swappable<T, T>::value)
    > {};

    template<class T, class U>
    struct is_nothrow_swappable : std::integral_constant<bool, 
        is_swappable<T, U>::value && noexcept(detail::swap_adl_tests::can_swap<T, U>(0))
    > {};
}

#endif
