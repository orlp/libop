#ifndef OP_UTILITY_H
#define OP_UTILITY_H

#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include "math.h"


namespace op {
    // Needed for declarations below.
    namespace detail {
        using swallow_param_pack_t = int[];
        template<class T, T N, int = (N == 0 || N == 1) ? N : 2> struct Seq;
    }

    // Expands and consumes parameter packs in left-to-right order.
    #define OP_SWALLOW_PARAM_PACK(...) \
        void(op::detail::swallow_param_pack_t{0, ((__VA_ARGS__), void(), 0)...})

    // make_array(...) is a helper function to create a std::array with a derived type and size.
    template<class... Types>
    constexpr std::array<typename std::common_type<Types...>::type, sizeof...(Types)>
    make_array(Types&&... t);

    // make_array<CT>(...) is a helper function to create a std::array with a derived size.
    template<class CT, class... Types>
    constexpr std::array<CT, sizeof...(Types)> make_array(Types&&... t);

    // to_array(c_array) is a helper function to create a std::array from a C array.
    template<class T, std::size_t N>
    constexpr std::array<T, N> to_array(T const (&a)[N]);

    // range(start, stop, step=1) returns an iterable object that iterates over {start, start +
    // step, start + 2*step, ...} as long as start + k*step < stop.
    class range;
    
    // A helper for tuple_visit, an object of this type gets passed into the callback function to
    // determine the default for the return type of the callback function.
    class tuple_visit_ret;

    // tuple_visit(tuple, index, func) returns func(std::get<index>(tuple)). std::out_of_range is
    // thrown when the index is out of bounds of the tuple.
    template<class Func, class Tuple,
             class Ret = decltype(std::declval<Func>()(std::declval<op::tuple_visit_ret>()))>
    Ret tuple_visit(Tuple&& tuple, std::size_t index, const Func& func);

    // visit_forward<T>()(arg) returns arg implicitly converted to T, throws a std::invalid_argument
    // if it can't be implicitly converted.
    template<class T> struct visit_forward;

    // Drop-in replacement for C++14 std::integer_sequence.
    template<class T, T... I>
    struct integer_sequence {
        using value_type = T;
        static constexpr std::size_t size() noexcept { return sizeof...(I); }
    };

    template<class T, T N> using make_integer_sequence = typename detail::Seq<T, N>::type;
    template<std::size_t... I> using index_sequence = integer_sequence<std::size_t, I...>;
    template<std::size_t N> using make_index_sequence = make_integer_sequence<std::size_t, N>;
    template<class... T> using index_sequence_for = make_index_sequence<sizeof...(T)>;
}



// Implementation.
namespace op {
    namespace detail {
        template<class T, std::size_t N, std::size_t ...I>
        inline constexpr std::array<T, N> array_unpacker(T (&a)[N], op::index_sequence<I...>) {
            return {{ a[I]... }};
        }


        // C++14's integer_sequence log(n) implementation.
        template<class, class> struct Combine;
        template<class T, T ...I, T ...J>
        struct Combine<integer_sequence<T, I...>, integer_sequence<T, J...>> {
            using type = integer_sequence<T, I..., (sizeof...(I) + J)...>;
        };
        
        template<class T, T N, int> // Default argument defined at start of file.
        struct Seq {
            static_assert(N >= 0, "N must be nonnegative in make_integer_sequence<T, N>");
            using type = typename Combine<typename Seq<T,     N/2>::type,
                                          typename Seq<T, N - N/2>::type>::type;
        };
        
        template<class T, T N> struct Seq<T, N, 0> { using type = integer_sequence<T>; };
        template<class T, T N> struct Seq<T, N, 1> { using type = integer_sequence<T, 0>; };


        template<class Ret, class Tuple, class Func>
        inline Ret tuple_visit_helper(Tuple&&, std::size_t, const Func&, op::index_sequence<>) {
            throw std::out_of_range("tuple index is out of bounds");
        }

        template<class Ret, class Tuple, class Func, std::size_t N, std::size_t... Tail>
        inline Ret tuple_visit_helper(Tuple&& tuple, std::size_t index, const Func& func,
                                      op::index_sequence<N, Tail...>) {
            if (index == N) return func(std::get<N>(tuple));

            return tuple_visit_helper<Ret>(std::forward<Tuple>(tuple), index, func,
                                           op::index_sequence<Tail...>());
        }
    }


    template<class... Types>
    inline constexpr std::array<typename std::common_type<Types...>::type, sizeof...(Types)>
    make_array(Types&&... t) {
        return {{ std::forward<Types>(t)... }};
    }


    template<class CT, class... Types>
    inline constexpr std::array<CT, sizeof...(Types)> make_array(Types&&... t) {
        return {{ static_cast<CT>(std::forward<Types>(t))... }};
    }


    template<class T, std::size_t N>
    inline constexpr std::array<T, N> to_array(T const (&a)[N]) {
        return detail::array_unpacker(a, op::make_index_sequence<N>());
    }


    class range {
    public:
        class iterator : public std::iterator<std::bidirectional_iterator_tag, int64_t> {
        public:
            int64_t operator*() const { return cur_; }
            iterator& operator++() { cur_ += step_; return *this; }
            iterator operator++(int) { iterator tmp(*this); cur_ += step_; return tmp; }
            iterator& operator--() { cur_ -= step_; return *this; }
            iterator operator--(int) { iterator tmp(*this); cur_ -= step_; return tmp; }

            bool operator==(const iterator& other) { return cur_ == other.cur_; }
            bool operator!=(const iterator& other) { return cur_ != other.cur_; }

        private:
            friend class range;

            iterator(int64_t cur, int64_t step) : cur_(cur), step_(step) { }

            int64_t cur_;
            int64_t step_;
        };

        template<class T>
        range(T stop) : start_(0), stop_(conv_int64(stop)), step_(1) { }

        template<class T, class U>
        range(T start, U stop) : start_(conv_int64(start)), stop_(conv_int64(stop)), step_(1) {
            if (stop_ < start_) start_ = stop_;
        }
        
        template<class T, class U, class V>
        range(T start, U stop, V step)
        : start_(conv_int64(start)), stop_(conv_int64(stop)), step_(conv_int64(step)) {
            if (step_ == 0) start_ = stop_;
            else if (step_ > 0) {
                if (stop_ < start_) start_ = stop_;

                auto size = stop_ - start_;
                stop_ = start_ + (1 + (size - 1) / step_) * step_;
            } else {
                if (start_ < stop_) start_ = stop_;

                auto size = start_ - stop_;
                stop_ = start_ - (1 + (size - 1) / -step_) * -step_;
            }
        }

        iterator begin() const { return {start_, step_}; }
        iterator end() const { return {stop_, step_}; }

    private:
        template <typename T>
        typename std::enable_if<
            op::safe_greater(std::numeric_limits<T>::max(), std::numeric_limits<int64_t>::max()),
        int64_t>::type conv_int64(T x) {
            if (x > static_cast<T>(std::numeric_limits<int64_t>::max())) {
                throw std::range_error("range() parameter too large, can't fit in int64_t");
            }

            return x;
        }

        template <typename T>
        typename std::enable_if<
            !op::safe_greater(std::numeric_limits<T>::max(), std::numeric_limits<int64_t>::max()),
        int64_t>::type conv_int64(T x) {
            return x;
        }

        int64_t start_;
        int64_t stop_;
        int64_t step_;
    };


    template<class Func, class Tuple, class Ret>
    inline Ret tuple_visit(Tuple&& tuple, std::size_t index, const Func& func) {
        return detail::tuple_visit_helper<Ret>(std::forward<Tuple>(tuple), index, func,
               op::make_index_sequence<std::tuple_size<typename std::decay<Tuple>::type>::value>());
    }


    template<class T>
    struct visit_forward {
        template<class U>
        typename std::enable_if<std::is_convertible<U, T>::value, T>::type
        operator()(U&& u) const { return std::forward<U>(u); }


        template<class U>
        typename std::enable_if<!std::is_convertible<U, T>::value, T>::type
        operator()(U&& u) const {
            throw std::invalid_argument("visit_forward called with unforwardable type");
        }
    };
}

#endif
