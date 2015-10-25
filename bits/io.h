#ifndef OP_IO_H
#define OP_IO_H

#include <cstddef>
#include <ios>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "utility.h"

/*
    Formatting mini-language specification
    --------------------------------------

    Most below formatting functions take a format string and an arbitrary amount of arguments.
    The format string gets streamed character by character into the output stream, until a
    replacement field is encountered. A replacement field is enclosed by curly braces. This
    means you'll have to escape literal curly braces by doubling them: {{ and }}.

    Formatting is stateless and totally independent of the current options/flags of the output
    stream. All options/flags of the stream are saved before formatting takes place and are
    restored to their original values when done. Neither side can affect the other.

    When a replacement field is encountered it gets parsed using the following grammar:

        replacement_field ::= "{" [index] [":" format_spec] "}"
        index             ::= <integer>
        format_spec       ::= [fill][align][width]["." precision][flag]*
        fill              ::= "'" <any character>
        align             ::= "<" | ">" | "^"
        width             ::= <integer> | "*" <integer>
        precision         ::= <integer> | "*" <integer>
        flag              ::= "b" | "e" | "f" | "o" | "p" | "t" | "u" | "x" | "+"

    Every replacement field takes an argument and streams it into the output stream using the
    requested options/flags. If no options are specified it will use the default options of
    std::ostream. Options/flags do not carry over between replacement fields.

    First the index of the argument is determined. You can specify this explicitly, for example
    "{0}" to select the first argument, or omit it and it will automatically be derived as one
    plus the last used index (0 for the first use). So the following two format strings are
    identical:

        "{0} {1} {0} {1}"
                  ^ You can re-use arguments for multiple replacement fields.

        "{} {} {0} {}"

    Then comes an optional format spec, starting with a colon (:).
    
    If a fill is specified (starting with an apostrophe (') followed by the fill character) it
    will be passed to out.fill(fill_char). The default fill is a space character.

    If an align is specified the following manipulator will be streamed into the output stream:

        < = std::left (this is the default)
        ^ = std::internal
        > = std::right

    Then it's possible to specify a width as an integer, passed to out.width(width). However, if
    the width parameter starts with an asterisk (*) the actual width is derived differently. You
    can specify an index immediately after the asterisk, where the formatting function will look
    for the width as an integer in the passed arguments. If the index is omitted it will
    automatically be derived as one plus the last used index. Note that the index for width
    counts as "the last used index" for future indexing operations. For example the three
    following format calls are all identical:

        format("{0:'_10} {1}\n", 12345, 42)
        format("{0:'_*1} {2}\n", 12345, 10, 42)
        format("{:'_*} {}\n",    12345, 10, 42)

        All result in "_____12345".
    
    If the width starts with a 0 and there is no other fill character specified, the fill character
    will be set to '0'. After the width the precision is parsed in similar fashion. The precision
    always starts with a period (.) and gets passed to out.precision(precision). The default
    precision is 6.

    Finally you can specify zero or more flags. Each flag will result in a manipulator being
    streamed into the output stream before the argument. The meaning of the flag codes is as
    following:

        t -> std::boolalpha
        f -> std::fixed
        x -> std::hex
        o -> std::oct
        e -> std::scientific
        b -> std::showbase
        p -> std::showpoint
        + -> std::showpos
        u -> std::uppercase
*/


namespace op {
    // Acts like std::cout << arg1 << " " << arg2 << " " << ... << " " << argn << "\n";.
    template<class... Args> void print(Args&&... args);
    
    // Acts like out << arg1 << " " << arg2 << " " << ... << " " << argn << "\n";.
    template<class Char, class Trait, class... Args>
    void fprint(std::basic_ostream<Char, Trait>& out, Args&&... args);

    // Uses the above formatting language to format args using format string str into out.
    template<class Char, class Trait, class... Args>
    void fprintf(std::basic_ostream<Char, Trait>& out, const Char* str, Args&&... args);

    // Uses fprintf on std::cout.
    template<class... Args>
    void printf(const char* str, Args&&... args);
    
    // Uses fprintf on a std::basic_ostringstream<Char, Trait, Alloc> and returns the result.
    template<class Char, class Trait = std::char_traits<Char>, class Alloc = std::allocator<Char>,
             class... Args>
    std::basic_string<Char, Trait, Alloc> format(const Char* str, Args&&... args);
}



// Implementation
namespace op {
    template<class Char, class Trait>
    inline void fprint(std::basic_ostream<Char, Trait>& out) { out << out.widen('\n'); }
    
    template<class Char, class Trait, class Arg, class... Args>
    inline void fprint(std::basic_ostream<Char, Trait>& out, Arg&& arg, Args&&... args) {
        out << std::forward<Arg>(arg);
        OP_SWALLOW_PARAM_PACK(out << out.widen(' ') << std::forward<Args>(args));
        out << out.widen('\n');
    }
    
    template<class... Args>
    inline void print(Args&&... args) { op::fprint(std::cout, std::forward<Args>(args)...); }


    namespace detail {
        // Changes n to contain the integer parsed from str if an integer is found at the start of
        // str, otherwise n is unchanged. Changes str to one past the last digit. Returns whether
        // an integer was parsed or not.
        template<class T, class Char>
        inline bool parse_integer(T& n, const Char*& str, const Char zero) {
            if (*str >= zero && *str <= zero + 9) {
                n = 0;
                do n = 10*n + (*str++ - zero);
                while (*str >= zero && *str <= zero + 9);

                return true;
            }

            return false;
        }


        // Used to output tuple elements in the formatting code with tuple_visit.
        template<class Char, class Trait>
        class StreamerVisit {
        public:
            StreamerVisit(std::basic_ostream<Char, Trait>& out) : out(out) { }
            template<typename T> void operator()(T&& t) const { out << std::forward<T>(t); }

        private:
            std::basic_ostream<Char, Trait>& out;
        };


        // Used to save, modify and restore stream state formatting parameters.
        template<class Char, class Trait>
        struct StreamState {
            using fmtflags = std::ios_base::fmtflags;

            std::streamsize width;
            std::streamsize precision;
            Char fill;
            fmtflags flags;

            StreamState(const std::basic_ostream<Char, Trait>& s)
            : width(s.width()), precision(s.precision()), fill(s.fill()), flags(s.flags()) { }

            static StreamState default_for_stream(const std::basic_ostream<Char, Trait>& s) {
                StreamState default_state;
                default_state.width = 0;
                default_state.precision = 6;
                default_state.fill = s.widen(' ');
                default_state.flags = s.dec | s.skipws;
                return default_state;
            }

            fmtflags setf(fmtflags new_flags) {
                auto old = flags;
                flags |= new_flags;
                return old;
            }

            fmtflags setf(fmtflags new_flags, fmtflags mask) {
                auto old = flags;
                flags = (flags & ~mask) | (new_flags & mask);
                return old;
            }

            void apply(std::basic_ostream<Char, Trait>& s) {
                s.width(width);
                s.precision(precision);
                s.fill(fill);
                s.flags(flags);
            }

        private:
            StreamState() = default;
        };
    }


    template<class Char, class Trait, class... Args>
    inline void fprintf(std::basic_ostream<Char, Trait>& out, const Char* str, Args&&... args) {
        // Re-use widen.
        const Char zero = out.widen('0');
        const Char open_bracket = out.widen('{');
        const Char close_bracket = out.widen('}');
        
        // Put arguments in tuple for easy extraction.
        const auto& args_tuple = std::make_tuple(std::forward<Args>(args)...);
        
        // Save original stream state.
        detail::StreamState<Char, Trait> original_state(out);

        // This auto_index will always be the last used index + 1. Whenever an index is omitted in
        // the format string auto_index will be used. 
        std::size_t auto_index = 0;

        while (*str) {
            // Output character if it's not special in the format string.
            if (!Trait::eq(*str, open_bracket)) {
                // Closing brackets should be escaped to balance escaped opening brackets.
                if (Trait::eq(*str, close_bracket) && !Trait::eq(*++str, close_bracket)) {
                    throw std::runtime_error("single '}' encountered in format string");
                }

                out << *str++; continue;
            }

            // We know *str == '{' now, increment and check for escape.
            if (Trait::eq(*++str, open_bracket)) {
                out << *str++; continue;
            }

            // Start of format parameter.
            auto fmt_params = detail::StreamState<Char, Trait>::default_for_stream(out);

            // Parse index.
            std::size_t arg_index = auto_index;
            detail::parse_integer(arg_index, str, zero);
            auto_index = arg_index + 1;

            // Parse format spec.
            if (Trait::eq(*str, out.widen(':'))) {
                ++str;

                // Check for fill.
                bool fill_specified = false;
                if (Trait::eq(*str, out.widen('\'')) && *(str + 1)) {
                    ++str;
                    fmt_params.fill = *str++;
                    fill_specified = true;
                }

                // Check for align.
                if (Trait::eq(*str, out.widen('<'))) {
                    ++str;
                    out.setf(out.left, out.adjustfield);
                } else if (Trait::eq(*str, out.widen('>'))) {
                    ++str;
                    out.setf(out.right, out.adjustfield);
                } else if (Trait::eq(*str, out.widen('^'))) {
                    ++str;
                    out.setf(out.internal, out.adjustfield);
                }

                // Check for width.
                if (Trait::eq(*str, out.widen('*'))) {
                    ++str;

                    std::size_t width_index = auto_index;
                    detail::parse_integer(width_index, str, zero);
                    auto_index = width_index + 1;

                    fmt_params.width = op::tuple_visit(args_tuple, width_index,
                                                       op::visit_forward<std::streamsize>());
                } else {
                    if (!fill_specified && Trait::eq(*str, zero)) fmt_params.fill = *str++;
                    detail::parse_integer(fmt_params.width, str, zero);
                }

                // Check for precision.
                if (Trait::eq(*str, out.widen('.'))) {
                    ++str;
  
                    if (Trait::eq(*str, out.widen('*'))) {
                        ++str;

                        std::size_t precision_index = auto_index;
                        detail::parse_integer(precision_index, str, zero);
                        auto_index = precision_index + 1;

                        fmt_params.precision =
                            op::tuple_visit(args_tuple, precision_index,
                                            op::visit_forward<std::streamsize>());
                    } else {
                        if (!detail::parse_integer(fmt_params.precision, str, zero)) {
                            throw std::runtime_error("expected precision after .");
                        }
                    }
                }
                
                // Check for flags.
                while (*str) {
                         if (Trait::eq(*str, out.widen('t'))) fmt_params.setf(out.boolalpha);
                    else if (Trait::eq(*str, out.widen('f'))) fmt_params.setf(out.fixed, out.floatfield);
                    else if (Trait::eq(*str, out.widen('x'))) fmt_params.setf(out.hex, out.basefield);
                    else if (Trait::eq(*str, out.widen('o'))) fmt_params.setf(out.oct, out.basefield);
                    else if (Trait::eq(*str, out.widen('e'))) fmt_params.setf(out.scientific, out.floatfield);
                    else if (Trait::eq(*str, out.widen('b'))) fmt_params.setf(out.showbase);
                    else if (Trait::eq(*str, out.widen('p'))) fmt_params.setf(out.showpoint);
                    else if (Trait::eq(*str, out.widen('+'))) fmt_params.setf(out.showpos);
                    else if (Trait::eq(*str, out.widen('u'))) fmt_params.setf(out.uppercase);
                    else break; // Unknown flag.
                    ++str;
                }
            }

            if (Trait::eq(*str, close_bracket)) {
                ++str;

                // Stream argument into out.
                fmt_params.apply(out);
                op::tuple_visit(args_tuple, arg_index, detail::StreamerVisit<Char, Trait>(out));
            } else {
                throw std::runtime_error("expected end of format string");
            }
        }

        // Restore original stream state.
        original_state.apply(out);
    }


    template<class Char, class Trait, class Alloc, class... Args>
    inline std::basic_string<Char, Trait, Alloc> format(const Char* str, Args&&... args) {
        std::basic_ostringstream<Char, Trait, Alloc> stream;
        op::fprintf(stream, str, std::forward<Args>(args)...);
        return stream.str();
    }


    template<class... Args>
    inline void printf(const char* str, Args&&... args) {
        op::fprintf(std::cout, str, std::forward<Args>(args)...);
    }
}

#endif
