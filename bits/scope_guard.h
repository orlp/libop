#ifndef OP_SCOPE_GUARD_H
#define OP_SCOPE_GUARD_H

#include <type_traits>
#include <utility>

#include "preprocessor.h"


/*
    Usage:

        OP_SCOPE_EXIT {
            // ...
        };

    Function body may not throw. Function body is at scope exit, regardless of exceptions thrown.
*/



// Implementation
#define OP_SCOPE_EXIT auto OP_ANONYMOUS_VAR(op_scope_exit_holder) = \
    ::op::detail::ScopeExitAcceptor() + [&]() noexcept 


namespace op {
    namespace detail {
        template<class F>
        struct ScopeExit {
            ScopeExit(const ScopeExit&) = delete;
            ScopeExit(ScopeExit&&) = default;
            ~ScopeExit() noexcept { fn(); }
            F fn;
        };

        struct ScopeExitAcceptor { };

        template<class F>
        ScopeExit<typename std::decay<F>::type> operator+(ScopeExitAcceptor, F&& fn) {
            return {std::forward<F>(fn)};
        }
    }
}

#endif
