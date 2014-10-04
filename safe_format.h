#ifndef SAFE_FORMAT_H
#define SAFE_FORMAT_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace safe_format {
    namespace detail {
        typedef std::string::const_iterator str_citer;

        // parses an integer from the options
        str_citer parse_int(str_citer options_begin, str_citer options_end, int& integer) {
            if (options_begin != options_end && std::isdigit(*options_begin)) {
                char* end;
                
                integer = std::strtol(&*options_begin, &end, 10);
                if (integer == std::numeric_limits<int>::max()) {
                    throw std::runtime_error("format string contains integer that is too large");
                }
                    
                options_begin += end - &*options_begin; 
            }
                
            return options_begin;
        }

        // parses options for a fill character and alignment and returns an iterator to one past the last parsed character
        str_citer parse_padding(str_citer options_begin, str_citer options_end, char& fill, char& align) {
            int len = std::distance(options_begin, options_end);
            if (len >= 2 && std::string("<>^=").find(options_begin[1]) != std::string::npos) {
                fill = options_begin[0];
                align = options_begin[1];            
                return options_begin + 2; 
            } else if (len >= 1 && std::string("<>^=").find(options_begin[0]) != std::string::npos) {
                align = options_begin[0];
                return options_begin + 1;
            }

            return options_begin;
        }

        // parses the sign and alternate options, and detects if leading zeroes should be used
        str_citer parse_numeric(str_citer options_begin, str_citer options_end, char& fill, char& align, char& sign, bool& alternate) {
            if (options_begin != options_end && (*options_begin == '+' || *options_begin == '-' || *options_begin == ' ')) {
                sign = *options_begin;
                ++options_begin;
            }

            if (options_begin != options_end && *options_begin == '#') {
                alternate = true;
                ++options_begin;
            }

            // look ahead to see if the width starts with a zero
            if (options_begin != options_end && *options_begin == '0') {
                fill = '0';
                align = '=';
            }

            return options_begin;
        }

        // parses the minimal width and precision from the options, if present
        str_citer parse_width_precision(str_citer options_begin, str_citer options_end, int& min_width, int& precision) {
            options_begin = parse_int(options_begin, options_end, min_width);
            
            if (options_begin != options_end && *options_begin == '.') {
                ++options_begin;
                options_begin = parse_int(options_begin, options_end, precision);
            }

            return options_begin;
        }
        
        // pads the given value properly using the given fill, alignment and minimal width
        std::string pad(const std::string& value, char fill, char align, int min_width) {
            int padding = min_width - value.length();
            if (padding > 0) {
                if (align == '<') {
                    return value + std::string(padding, fill);
                } else if (align == '^') {
                    return std::string(padding / 2, fill) + value + std::string(padding - padding / 2, fill);
                } else if (align == '>') {
                    return std::string(padding, fill) + value;
                } else if (align == '=') {
                    if (!value.empty() && (value[0] == '+' || value[0] == '-')) {
                        return value[0] + std::string(padding, fill) + value;
                    }

                    return std::string(padding, fill) + value;
                }
            }

            return value;
        }

        // a general fallback implementation using stringstreams
        template<class T>
        typename std::enable_if<!std::is_integral<T>::value, std::string>::type safe_format_value(str_citer options_begin, str_citer options_end, const T& arg) {
            char fill = ' ';
            char align = '<';
            int min_width = 0;
            int precision = -1;

            options_begin = parse_padding(options_begin, options_end, fill, align);
            options_begin = parse_width_precision(options_begin, options_end, min_width, precision);

            if (precision != -1) {
                throw std::runtime_error("precision not valid for fallback implementation");
            }

            if (options_begin != options_end) {
                throw std::runtime_error("format string too long, expected end");
            }

            std::stringstream s;
            s << arg;
            return pad(s.str(), fill, align, min_width);
        }

        // an implementation for integral types
        template<class T>
        typename std::enable_if<std::is_integral<T>::value, std::string>::type safe_format_value(str_citer options_begin, str_citer options_end, const T& arg) {
            char fill = ' ';
            char align = '<';
            char sign = '-';
            bool alternate = false;
            int min_width = 0;
            int precision = -1;

            options_begin = parse_padding(options_begin, options_end, fill, align);
            options_begin = parse_numeric(options_begin, options_end, fill, align, sign, alternate);
            options_begin = parse_width_precision(options_begin, options_end, min_width, precision);

            std::string result;
            const char digits[] = "0123456789abcdef";
            const char DIGITS[] = "0123456789ABCDEF";

            T conv(arg);
            while (conv > 0) {
                result += digits[conv % 10];
                conv /= 10;
            }

            return "int";
        }

        // if we reach this function through recursion we've ran out of arguments, so the index was out of range
        std::string format_value_unpacker(int idx, str_citer options_begin, str_citer options_end) {
            throw std::runtime_error("format argument index out of range");
        }

        // recursively find the correct argument to format
        template<class T, class... Args>
        std::string format_value_unpacker(int idx, str_citer options_begin, str_citer options_end, const T& arg, Args... args) {
            if (idx == 0) return safe_format_value(options_begin, options_end, arg);
            return format_value_unpacker(--idx, options_begin, options_end, args...);
        }
    }

    // The format function, splits the format string up for each individual argument and then
    // dispatches them to individual formatters with the correct options.
    template<class... Args>
    std::string format(const std::string& fmt, Args... args) {
        std::string result;
        int arg_index = -1;
        auto fmt_it = fmt.cbegin();

        while (fmt_it != fmt.end()) {
            // check for opening bracket
            if (*fmt_it == '{') {
                if (++fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");
                
                // double bracket escapes
                if (*fmt_it == '{') {
                    result += '{';
                    ++fmt_it;

                // it's an actual field
                } else {
                    // is an index specified or use an automatic index?
                    if (std::isdigit(*fmt_it)) {
                        fmt_it = detail::parse_int(fmt_it, fmt.end(), arg_index);
                        if (fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");
                    } else {
                        ++arg_index;
                    }
                    
                    // save the current index, this might get overwritten by nested brackets in the options
                    int unnested_arg_index = arg_index;

                    // are there options?
                    std::string options;
                    if (*fmt_it == ':') {
                        if (++fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");

                        // get options, handling nested brackets
                        while (*fmt_it != '}') {
                            // nested brackets
                            if (*fmt_it == '{') {
                                if (++fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");

                                if (std::isdigit(*fmt_it)) {
                                    fmt_it = detail::parse_int(fmt_it, fmt.end(), arg_index);
                                    if (fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");
                                } else {
                                    ++arg_index;
                                }

                                if (*fmt_it != '}') throw std::runtime_error("unmatched '{' in format");
                                if (++fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");
                                options += detail::format_value_unpacker(arg_index, "", args...);
                            } else {
                                options += *fmt_it;
                            }
                        
                            if (++fmt_it == fmt.end()) throw std::runtime_error("unmatched '{' in format");
                        }
                    }
        
                    result += detail::format_value_unpacker(unnested_arg_index, options, args...);
                }

            // closing bracket, test for escape or error
            } else if (*fmt_it == '}') {
                if (++fmt_it == fmt.end() || *fmt_it != '}') {
                    throw std::runtime_error("single '}' encountered in format");
                }

                result += '}';
                ++fmt_it;

            // any other character goes directly into the result
            } else {
                result += *fmt_it++;
            }
        }

        return result;
    }
}

using safe_format::format;

#endif
