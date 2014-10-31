#ifndef OP_STRING_H
#define OP_STRING_H

#include <cctype>
#include <cstddef>
#include <ostream>
#include <set>
#include <sstream>
#include <string>


namespace op {
    // Returns a string that is equal to every element in the sequence begin, end joined by sep.
    // Uses operator<<(std::ostream, T) to turn the elements into strings if they aren't already.
    template<class InIter>
    std::string join(InIter begin, InIter end, const std::string& sep = " ");

    // Inspects str and calls *out++ = token for each contiguous sequence of characters found in
    // str that does not contain any whitespace (determined by std::isspace). So " 1  2" will
    // result in *out++ = "1"; *out++ = "2". Returns the iterator one past the last element
    // written.
    template<class OutIter>
    OutIter split(const std::string& str, OutIter out);

    // Splits str up in tokens seperated by sep, even if the token is empty.
    // split("1::2", ":", out) will result in *out++ = "1"; *out++ = ""; *out++ = "2";
    // Returns the iterator one past the last element written.
    template<class OutIter>
    OutIter split(const std::string& str, const std::string& sep, OutIter out);

    // Removes characters from the left of str as long as they're found in chars.
    // Example: lstrip("xxIx", "x") -> "Ix".
    std::string lstrip(const std::string& str, const std::string& chars = " \f\n\r\t\v");

    // Removes characters from the right of str as long as they're found in chars.
    // Example: rstrip("xxIx", "x") -> "xxI".
    std::string rstrip(const std::string& str, const std::string& chars = " \f\n\r\t\v");

    // Removes characters from the both sides of str as long as they're found in chars.
    // Example: strip("xxIx", "x") -> "I".
    std::string strip(const std::string& str, const std::string& chars = " \f\n\r\t\v");
    
    // Does str start with start?
    bool startswith(const std::string& str, const std::string& start);

    // Does str end with end?
    bool endswith(const std::string& str, const std::string& end);
}



// Implementation.
namespace op {
    template<class InIter>
    inline std::string join(InIter begin, InIter end, const std::string& sep) {
        std::stringstream s;

        if (begin != end) s << *begin;
        for (++begin; begin != end; ++begin) s << sep << *begin;
        
        return s.str();
    }


    template<class OutIter>
    inline OutIter split(const std::string& str, OutIter out) {
        std::string token;

        for (auto c : str) {
            if (std::isspace(c)) {
                if (token.size()) *out++ = token;
                token = "";
            } else token += c;
        }

        if (token.size()) *out++ = token;

        return out;
    }


    template<class OutIter>
    inline OutIter split(const std::string& str, const std::string& sep, OutIter out) {
        std::set<char> sep_key(sep.begin(), sep.end());
        std::string token;

        for (auto c : str) {
            if (sep_key.count(c)) {
                *out++ = token;
                token = "";
            } else token += c;
        }

        *out++ = token;

        return out;
    }


    inline std::string lstrip(const std::string& str, const std::string& chars) {
        std::size_t lpos = str.find_first_not_of(chars);

        if (lpos == std::string::npos) return "";
        return str.substr(lpos);
    }


    inline std::string rstrip(const std::string& str, const std::string& chars) {
        std::size_t rpos = str.find_last_not_of(chars);
        
        if (rpos == std::string::npos) return "";
        return str.substr(0, rpos + 1);
    }


    inline std::string strip(const std::string& str, const std::string& chars) {
        std::size_t lpos = str.find_first_not_of(chars);

        if (lpos == std::string::npos) return "";
        return str.substr(lpos, str.find_last_not_of(chars) - lpos + 1);
    }
    

    inline bool startswith(const std::string& str, const std::string& start) {
        return str.compare(0, start.length(), start) == 0;
    }


    inline bool endswith(const std::string& str, const std::string& end) {
        return end.length() <= str.length() && str.compare(str.length() - end.length(), end.length(), end) == 0;
    }
}

#endif
