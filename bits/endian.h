/*
    Copyright Redshift Software, Inc. 2013

    Boost Software License - Version 1.0 - August 17th, 2003

    Permission is hereby granted, free of charge, to any person or organization
    obtaining a copy of the software and accompanying documentation covered by
    this license (the "Software") to use, reproduce, display, distribute,
    execute, and transmit the Software, and to prepare derivative works of the
    Software, and to permit third-parties to whom the Software is furnished to
    do so, all subject to the following:

    The copyright notices in the Software and this entire statement, including
    the above license grant, this restriction and the following disclaimer,
    must be included in all copies of the Software, in whole or in part, and
    all derivative works of the Software, unless such copies or derivative
    works are solely in the form of machine-executable object code generated by
    a source language processor.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
    SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
    FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

#ifndef OP_ENDIAN_H
#define OP_ENDIAN_H

/* GNU libc provides a header defining __BYTE_ORDER, or _BYTE_ORDER.
   And some OS's provide some for of endian header also. */
#if defined(__GLIBC__) || defined(__GNU_LIBRARY__)
    #include <endian.h>
#else
    #if defined(macintosh) || defined(Macintosh) || \
       (defined(__APPLE__) && defined(__MACH__))
        #include <machine/endian.h>
    #else
        #if defined(BSD) || defined(_SYSTYPE_BSD)
            #if defined(__OpenBSD__)
                #include <machine/endian.h>
            #else
                #include <sys/endian.h>
            #endif
        #endif
    #endif
#endif
#if defined(__BYTE_ORDER)
    #if (__BYTE_ORDER == __BIG_ENDIAN)
        #define OP_ENDIAN_BIG
    #endif
    #if (__BYTE_ORDER == __LITTLE_ENDIAN)
        #define OP_ENDIAN_LITTLE
    #endif
    #if (__BYTE_ORDER == __PDP_ENDIAN)
        #define OP_ENDIAN_PDP
    #endif
#endif
#if !defined(__BYTE_ORDER) && defined(_BYTE_ORDER)
    #if (_BYTE_ORDER == _BIG_ENDIAN)
        #define OP_ENDIAN_BIG
    #endif
    #if (_BYTE_ORDER == _LITTLE_ENDIAN)
        #define OP_ENDIAN_LITTLE
    #endif
    #if (_BYTE_ORDER == _PDP_ENDIAN)
        #define OP_ENDIAN_PDP
    #endif
#endif

/* Big-endian. */
#if !defined(OP_ENDIAN_BIG) && !defined(OP_ENDIAN_LITTLE) && \
    !defined(OP_ENDIAN_PDP)
    #if (defined(__BIG_ENDIAN__) && !defined(__LITTLE_ENDIAN__)) || \
         defined(__ARMEB__) || \
         defined(__THUMBEB__) || \
         defined(__AARCH64EB__) || \
         defined(_MIPSEB) || \
         defined(__MIPSEB) || \
         defined(__MIPSEB__)
        #define OP_ENDIAN_BIG
   #endif
#endif

/* Little-endian. */
#if !defined(OP_ENDIAN_BIG) && !defined(OP_ENDIAN_LITTLE) && \
    !defined(OP_ENDIAN_PDP)
    #if (defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)) || \
         defined(__ARMEL__) || \
         defined(__THUMBEL__) || \
         defined(__AARCH64EL__) || \
         defined(_MIPSEL) || \
         defined(__MIPSEL) || \
         defined(__MIPSEL__)
        #define OP_ENDIAN_LITTLE
   #endif
#endif

/* Some architectures are strictly one endianess (as opposed
   the current common bi-endianess). */
#if !defined(OP_ENDIAN_BIG) && !defined(OP_ENDIAN_LITTLE) && \
    !defined(OP_ENDIAN_PDP)
    #if defined(__m68k__) || defined(M68000) || \
        defined(__hppa__) || defined(__hppa) || defined(__HPPA__) || \
        defined(__370__) || defined(__THW_370__) || \
        defined(__s390__) || defined(__s390x__) || \
        defined(__SYSC_ZARCH__)
        #define OP_ENDIAN_BIG
    #endif
    #if defined(__amd64__) || defined(__amd64) || \
        defined(__x86_64__) || defined(__x86_64) || \
        defined(_M_X64)  || \
        defined(__ia64__) || defined(_IA64) || \
        defined(__IA64__) || defined(__ia64) || \
        defined(_M_IA64) || defined(__itanium__) || \
        defined(i386) || defined(__i386__) || \
        defined(__i486__) || defined(__i586__) || \
        defined(__i686__) || defined(__i386) || \
        defined(_M_IX86) || defined(_X86_) || \
        defined(__THW_INTEL__) || defined(__I86__) || \
        defined(__INTEL__) || \
        defined(__bfin__) || defined(__BFIN__) || \
        defined(bfin) || defined(BFIN)
        #define OP_ENDIAN_LITTLE
    #endif
#endif

/* Windows on ARM, if not otherwise detected/specified, is
   always little-endian. */
#if !defined(OP_ENDIAN_BIG) && !defined(OP_ENDIAN_LITTLE) && \
    !defined(OP_ENDIAN_PDP)
    #if defined(__arm__) || defined(__thumb__) || \
        defined(__TARGET_ARCH_ARM) || defined(__TARGET_ARCH_THUMB)
       #if defined(_WIN32) || defined(_WIN64) || \
           defined(__WIN32__) || defined(__TOS_WIN__) || \
           defined(__WINDOWS__)
            #define OP_ENDIAN_LITTLE
       #endif
   #endif
#endif

#if defined(OP_ENDIAN_LITTLE)
    #define OP_ENDIAN_BYTE_ORDER 1234
#elif defined(OP_ENDIAN_BIG)
    #define OP_ENDIAN_BYTE_ORDER 4321
#elif defined(OP_ENDIAN_PDP)
    #define OP_ENDIAN_BYTE_ORDER 2134
#else
    #error Could not detect endianness.
#endif

#include <algorithm>

namespace op {
    // TODO: optimized implementations for small sizes
    template<class T>
    T byteswap(T x) {
        std::reverse(reinterpret_cast<unsigned char*>(&x), reinterpret_cast<unsigned char*>(&x) + sizeof(x));
        return x;
    }

    template<class T>
    T htole(T x) {
        #if defined(OP_ENDIAN_LITTLE)
            return x;
        #elif defined(OP_ENDIAN_BIG)
            return byteswap(x);
        #else
            #error PDP endian not supported.
        #endif
    }

    template<class T>
    T htobe(T x) {
        #if defined(OP_ENDIAN_LITTLE)
            return byteswap(x);
        #elif defined(OP_ENDIAN_BIG)
            return x;
        #else
            #error PDP endian not supported.
        #endif
    }

    template<class T>
    T letoh(T x) {
        return htole(x);
    }

    template<class T>
    T betoh(T x) {
        return htobe(x);
    }
}

#endif

