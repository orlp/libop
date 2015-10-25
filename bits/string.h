#ifndef OP_STRING_H
#define OP_STRING_H

#include <cctype>
#include <cstddef>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>

#include "type_traits.h"


namespace op {
    // Used for declarations below.
    namespace detail {
        template<class T, class = void>
        struct underlying_char_type_impl { };
        
        template<class T>
        struct underlying_char_type_impl<T*,
            typename std::enable_if<op::is_one_of<typename std::decay<T>::type,
                char, signed char, unsigned char, wchar_t, char16_t, char32_t>{}>::type> {
            typedef typename std::decay<T>::type type;
        };
        
        template<class CharT, class Traits, class Allocator>
        struct underlying_char_type_impl<std::basic_string<CharT, Traits, Allocator>, void> {
            typedef typename std::decay<CharT>::type type;
        };
        
        template<class T> using underlying_char_type =
            typename underlying_char_type_impl<typename std::decay<T>::type>::type;
        template<class... T> using combinable_string_types =
            op::is_all_same<underlying_char_type<typename std::decay<T>::type>...>;
        template<class... T> using enable_if_combinable_string_types_t =
            typename std::enable_if<combinable_string_types<T...>{}>::type;
        template<class T> using as_basic_string =
            std::basic_string<underlying_char_type<typename std::decay<T>::type>>;
    }

    // Returns a string that is equal to every element in the sequence begin, end joined by sep.
    // Uses operator<<(std::ostream, T) to turn the elements into strings if they aren't already.
    template<class InIter, class Str = const char*,
             class = detail::enable_if_combinable_string_types_t<Str>>
    detail::as_basic_string<Str> join(InIter begin, InIter end, const Str& sep = " ");

    // Inspects str and calls *out++ = token for each contiguous sequence of characters found in
    // str that does not contain any whitespace (determined by std::isspace). So " 1  2" will
    // result in *out++ = "1"; *out++ = "2". Returns the iterator one past the last element
    // written.
    template<class OutIter, class Str,
             class = detail::enable_if_combinable_string_types_t<Str>>
    OutIter split(const Str& str, OutIter out);

    // Splits str up in tokens separated by sep, even if the token is empty.
    // split("::1::2", "::", out) will result in *out++ = ""; *out++ = "1"; *out++ = "2";
    // Returns the iterator one past the last element written.
    template<class OutIter, class Str1, class Str2,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    OutIter split(const Str1& str, const Str2& sep, OutIter out);
    
    // Splits str up in tokens separated by one of the characters in sep, even if the token is empty.
    // split("+1-2", "+-", out) will result in *out++ = ""; *out++ = "1"; *out++ = "2";
    // Returns the iterator one past the last element written.
    template<class OutIter, class Str1, class Str2,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    OutIter split_one_of(const Str1& str, const Str2& sep, OutIter out);

    // Removes characters from the left of str as long as they're found in chars.
    // Example: lstrip("xxIx", "x") -> "Ix".
    template<class Str1, class Str2 = const char*,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    detail::as_basic_string<Str1> lstrip(const Str1& str, const Str2& chars = " \f\n\r\t\v");

    // Removes characters from the right of str as long as they're found in chars.
    // Example: rstrip("xxIx", "x") -> "xxI".
    template<class Str1, class Str2 = const char*,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    detail::as_basic_string<Str1> rstrip(const Str1& str, const Str2& chars = " \f\n\r\t\v");

    // Removes characters from the both sides of str as long as they're found in chars.
    // Example: strip("xxIx", "x") -> "I".
    template<class Str1, class Str2 = const char*,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    detail::as_basic_string<Str1> strip(const Str1& str, const Str2& chars = " \f\n\r\t\v");
    
    // Does str start with start?
    template<class Str1, class Str2,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    bool startswith(const Str1& str, const Str2& start);

    // Does str end with end?
    template<class Str1, class Str2,
             class = detail::enable_if_combinable_string_types_t<Str1, Str2>>
    bool endswith(const Str1& str, const Str2& end);
}



// Implementation.
namespace op {
    namespace detail {
        template<class CharT>
        inline std::basic_string<CharT, std::char_traits<CharT>, std::allocator<CharT>>
        to_str(const CharT* p)
            { return p; }
        
        template<class CharT, class Traits, class Allocator>
        inline const std::basic_string<CharT, Traits, Allocator>&
        to_str(const std::basic_string<CharT, Traits, Allocator>& s)
            { return s; }
    }

    template<class InIter, class Str, class>
    inline detail::as_basic_string<Str> join(InIter begin, InIter end, const Str& sep) {
        if (begin == end) return {};

        std::basic_stringstream<detail::underlying_char_type<Str>> s;
        s << *begin;
        for (++begin; begin != end; ++begin) s << sep << *begin;
        
        return s.str();
    }


    template<class OutIter, class Str, class>
    inline OutIter split(const Str& str, OutIter out) {
        detail::as_basic_string<Str> token;
        auto p = &str[0];

        while (*p) {
            auto c = *p++;
            if (std::isspace(c)) {
                if (token.size()) *out++ = token;
                token.clear();
            } else token.push_back(c);
        }

        if (token.size()) *out++ = token;

        return out;
    }


    template<class OutIter, class Str1, class Str2, class>
    inline OutIter split(const Str1& str_, const Str2& sep_, OutIter out) {
        const auto& str = detail::to_str(str_);
        const auto& sep = detail::to_str(sep_);
        std::size_t start = 0;
        std::size_t end = str.find(sep);

        while (end != str.npos) {
            auto token = str.substr(start, end - start);
            if (token.size()) *out++ = token;
            start = end + sep.length();
            end = str.find(sep, start);
        }

        if (start != str.size()) *out++ = str.substr(start, str.size());

        return out;
    }


    template<class OutIter, class Str1, class Str2, class>
    inline OutIter split_one_of(const Str1& str, const Str2& sep, OutIter out) {
        std::set<detail::underlying_char_type<Str1>> sep_key(sep.begin(), sep.end());
        detail::as_basic_string<Str1> token;

        auto p = &str[0];
        while (*p) {
            auto c = *p++;
            if (sep_key.count(c)) {
                *out++ = token;
                token.clear();
            } else token.push_back(c);
        }

        *out++ = token;

        return out;
    }


    template<class Str1, class Str2, class>
    inline detail::as_basic_string<Str1> lstrip(const Str1& str_, const Str2& chars) {
        const auto& str = detail::to_str(str_);
        std::size_t lpos = str.find_first_not_of(chars);
        if (lpos == str.npos) return {};
        return str.substr(lpos);
    }


    template<class Str1, class Str2, class>
    inline detail::as_basic_string<Str1> rstrip(const Str1& str_, const Str2& chars) {
        const auto& str = detail::to_str(str_);
        std::size_t rpos = str.find_last_not_of(chars);
        if (rpos == str.npos) return {};
        return str.substr(0, rpos + 1);
    }


    template<class Str1, class Str2, class>
    inline detail::as_basic_string<Str1> strip(const Str1& str_, const Str2& chars) {
        const auto& str = detail::to_str(str_);
        std::size_t lpos = str.find_first_not_of(chars);
        if (lpos == str.npos) return {};
        return str.substr(lpos, str.find_last_not_of(chars) - lpos + 1);
    }
    

    template<class Str1, class Str2, class>
    inline bool startswith(const Str1& str_, const Str2& start_) {
        const auto& str = detail::to_str(str_);
        const auto& start = detail::to_str(start_);
        return str.compare(0, start.length(), start) == 0;
    }


    template<class Str1, class Str2, class>
    inline bool endswith(const Str1& str_, const Str2& end_) {
        const auto& str = detail::to_str(str_);
        const auto& end = detail::to_str(end_);
        return end.length() <= str.length() &&
               str.compare(str.length() - end.length(), end.length(), end) == 0;
    }
}

#endif
