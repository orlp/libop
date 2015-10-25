#ifndef OP_TYPE_TRAITS_H
#define OP_TYPE_TRAITS_H

#include <cstddef>
#include <type_traits>
#include <cmath>


namespace op {
    // Inherits from std::true_type if T and U are Swappable, std::false_type if otherwise.
    template<class T, class U = T>
    struct is_swappable;

    // Same as is_swappable, with the added condition that the swap must be noexcept.
    template<class T, class U = T>
    struct is_nothrow_swappable;

    // Common floating point type. Turns Integral into double, rejects non-arithmetic types and 
    // then takes the common type.
    template<class... T>
    struct common_floating_point_type;

    template<class... T>
    using common_floating_point_type_t = typename common_floating_point_type<T...>::type;

    // Checks if T is the same as one of Us.
    template<class T, class... Us>
    struct is_one_of;

    // Checks if all types in Ts are the same.
    template<class... Ts> struct is_all_same;
}



// Implementation.
namespace op {
    namespace detail {
        namespace swap_adl_tests {
            // If swap ADL finds this then it would call std::swap otherwise (same signature).
            struct tag {};

            template<class T> tag swap(T&, T&);
            template<class T, std::size_t N> tag swap(T (&a)[N], T (&b)[N]);

            // Helper functions to test if an unqualified swap is possible, and if it becomes std::swap.
            template<class, class> std::false_type can_swap(...) noexcept(false);
            template<class T, class U, class = decltype(swap(std::declval<T&>(), std::declval<U&>()))>
            std::true_type can_swap(int) noexcept(
                noexcept(swap(std::declval<T&>(), std::declval<U&>()))
            );

            template<class, class> std::false_type uses_std(...);
            template<class T, class U>
            std::is_same<decltype(swap(std::declval<T&>(), std::declval<U&>())), tag> uses_std(int);

            // Helper templates to determine if swapping is noexcept. The C++11/14 standards have a
            // broken noexcept specification for multidimensional arrays, so we use a template solution.
            template<class T>
            struct is_std_swap_noexcept : std::integral_constant<bool,
                std::is_nothrow_move_constructible<T>::value &&
                std::is_nothrow_move_assignable<T>::value
            > { };

            template<class T, std::size_t N>
            struct is_std_swap_noexcept<T[N]> : is_std_swap_noexcept<T> { };

            template<class T, class U>
            struct is_adl_swap_noexcept : std::integral_constant<bool, noexcept(can_swap<T, U>(0))> { };
        }

        template<class T, class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
        struct common_ftype_helper {
            typedef T type;
        };
        
        template<class T>
        struct common_ftype_helper<T, typename std::enable_if<std::is_integral<T>::value>::type> {
            typedef double type;
        };

        template<class T>
        using common_ftype_helper_t = typename common_ftype_helper<T>::type;
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
        is_swappable<T, U>::value && (
            (decltype(detail::swap_adl_tests::uses_std<T, U>(0))::value &&
                detail::swap_adl_tests::is_std_swap_noexcept<T>::value)
            ||
            (!decltype(detail::swap_adl_tests::uses_std<T, U>(0))::value &&
                detail::swap_adl_tests::is_adl_swap_noexcept<T, U>::value)
        )
    > {};

    template<class... T>
    struct common_floating_point_type {
        typedef typename std::common_type<detail::common_ftype_helper_t<T>...>::type type;
    };

    template<class T, class... Us>
    struct is_one_of : std::false_type {};

    template<class T, class U, class... Us>
    struct is_one_of<T, U, Us...> :
        std::integral_constant<bool, std::is_same<T, U>{} || is_one_of<T, Us...>{}> {};
    
    template<class... Ts> struct is_all_same : std::false_type {};
    template<> struct is_all_same<> : std::true_type {};
    template<class T> struct is_all_same<T> : std::true_type {};
    template<class T, class U, class... Us>
    struct is_all_same<T, U, Us...>
        : std::integral_constant<bool, std::is_same<T, U>{} && is_all_same<U, Us...>{}> {};
}

#endif



