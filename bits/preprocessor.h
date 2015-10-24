#ifndef OP_PREPROCESSOR_H
#define OP_PREPROCESSOR_H

#define OP_CONCAT(a, b) OP_CONCAT_EXPANDED(a, b)
#define OP_CONCAT_EXPANDED(a, b) a##b

#ifdef __COUNTER__
    #define OP_ANONYMOUS_VAR(prefix) OP_CONCAT(prefix, __COUNTER__)
#else
    #define OP_ANONYMOUS_VAR(prefix) OP_CONCAT(prefix, __LINE__)
#endif

#endif
