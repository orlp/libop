#ifndef OP_RANDOM_H
#define OP_RANDOM_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
    #include <windows.h>
    #include <wincrypt.h>
#else
    #include <fstream>
#endif

#include "exception.h"


namespace op {
    // Returns an iterator pointing to a random element from the given sequence.
    template<class Iter, class URNG>
    Iter random_choice(Iter first, Iter last, URNG&& g);

    // Selects k random samples without replacement from the given sequence. The order of the
    // samples in the output is undefined, but not necessarily random.
    template<class InIter, class OutIter, class URNG>
    void random_sample(InIter first, InIter last, OutIter out, int k, URNG&& g);

    // Fill buf with buflen random bytes. Slower on big-endian machines for consistency.
    template<class URNG>
    void random_bytes(unsigned char* buf, int buflen, URNG&& g);

    // randint<T>(g) returns a random integer of type T, with a uniform bit pattern.
    template<class T, class URNG>
    T randint(URNG&& g);

    // randint(min, max, g) returns a random integer N with min <= N <= max.
    template<class T, class URNG>
    T randint(T min, T max, URNG&& g);

    // rand(g) returns a random floating point number in the range [0.0, 1.0).
    template<class T = double, class URNG>
    T rand(URNG&& g);

    // rand(min, max, g) returns a random floating point number N such that min <= N <= max.
    template<class T, class URNG>
    T rand(T min, T max, URNG&& g);

    // Exception used by op::random_device.
    struct RandomDeviceError : public virtual op::BaseException {
        explicit RandomDeviceError(const std::string& msg) : op::BaseException(msg) { }
        explicit RandomDeviceError(const char* msg) : op::BaseException(msg) { }
    };

    // C++11 <random> compatible random_device, because many implementations are broken.
    // Implementation is cryptographically secure, or fails with a RandomDeviceError.
    class random_device {
    public:
        typedef uint32_t result_type;

        explicit random_device(const std::string& token =
            #ifdef _WIN32
                MS_DEF_PROV
            #else
                "/dev/urandom"
            #endif
        ) : token(token) { };

        random_device(const random_device&) = delete;
        random_device& operator=(random_device) = delete;

        result_type operator()();

        double entropy() const noexcept { return 8 * sizeof(result_type); }
        static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
        static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

    private:
        std::string token;
    };
}



// Implementation.
namespace op {
    namespace detail {
        template<class Iter, class URNG>
        inline Iter random_choice(Iter first, Iter last, URNG&& g, std::random_access_iterator_tag) {
            return first + op::randint<typename std::iterator_traits<Iter>::difference_type>(0, last - first - 1, g);
        }


        template<class Iter, class URNG>
        inline Iter random_choice(Iter first, Iter last, URNG&& g, std::forward_iterator_tag) {
            typedef typename std::iterator_traits<Iter>::difference_type diff_t;

            Iter result = first++;
            for (diff_t i = 1; first != last; ++i, ++first) {
                if (op::randint(0, i, g) == 0) result = first;
            }

            return result;
        }
    }


    template<class Iter, class URNG>
    inline Iter random_choice(Iter first, Iter last, URNG&& g) {
        if (first == last) throw std::out_of_range("empty sequence passed to random_choice");

        return detail::random_choice(first, last, g, typename std::iterator_traits<Iter>::iterator_category());
    }


    template<class InIter, class OutIter, class URNG>
    inline void random_sample(InIter first, InIter last, OutIter out, int k, URNG&& g) {
        typedef typename std::iterator_traits<InIter>::value_type T;
        typedef typename std::iterator_traits<InIter>::difference_type diff_t;

        std::vector<T> reservoir;
        reservoir.reserve(k);

        for (diff_t i = 0; i < k; ++i, ++first) {
            if (first == last) throw std::out_of_range("k larger than sequence size passed to random_sample");
            reservoir.push_back(*first);
        }

        for (diff_t i = k; first != last; ++i, ++first) {
            diff_t j = std::uniform_int_distribution<diff_t>(0, i)(g);
            if (j < k) reservoir[j] = *first;
        }

        for (auto& e : reservoir) *out++ = std::move(e);
    }


    template<class OutIter, class URNG>
    inline void random_bytes(OutIter out, size_t num, URNG&& g) {
        while (num--) *out++ = randint<unsigned char>(g);
    }


    template<class T, class URNG>
    inline T randint(URNG&& g) {
        return std::uniform_int_distribution<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max())(g);
    }


    template<class T, class URNG>
    inline T randint(T min, T max, URNG&& g) {
        return std::uniform_int_distribution<T>(min, max)(g);
    }


    template<class T, class URNG>
    inline T rand(URNG&& g) {
        return std::generate_canonical<T, std::numeric_limits<T>::digits>(g);
    }


    template<class T, class URNG>
    inline T rand(T min, T max, URNG&& g) {
        return std::uniform_real_distribution<T>(min, max)(g);
    }


    inline random_device::result_type op::random_device::operator()() {
        result_type result;

        #ifdef _WIN32
            HCRYPTPROV prov;

            if (!CryptAcquireContext(&prov, NULL, token.c_str(), PROV_RSA_FULL, CRYPT_VERIFYCONTEXT))  {
                throw op::RandomDeviceError("CryptAcquireContext failed");
            }

            if (!CryptGenRandom(prov, sizeof(result), reinterpret_cast<unsigned char*>(&result)))  {
                throw op::RandomDeviceError("CryptGenRandom failed");
            }

            CryptReleaseContext(prov, 0);
        #else
            std::ifstream device(token.c_str(), std::ios::in | std::ios::binary);

            if (!device.read(reinterpret_cast<char*>(&result), sizeof(result))) {
                throw op::RandomDeviceError("reading from the operating systems random device failed");
            }
        #endif

        return result;
    }
}

#endif
