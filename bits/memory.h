#ifndef OP_MEMORY_H
#define OP_MEMORY_H

#include <cassert>
#include <cstddef>
#include <map>
#include <new>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <vector>


namespace op {
    /*
        arena<N, T=char> holds enough stack storage for N objects of type T. The arena is not
        thread-safe, only use it on one thread. Its methods are:

            // Amount of bytes of storage that is allocated.
            std::size_t used() const noexcept; 

            // Size in bytes of the arena.
            static constexpr std::size_t size() noexcept;

            // Reclaims all allocated memory to get back to the initial state. Beware: arena only
            // handles uninitialized data - this does not call destructors.
            void reset() noexcept;  

            // Allocates and returns the address to n bytes from the arena, or returns nullptr if
            // there is not enough space left to allocate. The returned memory has the alignment of
            // std::max_align_t.
            void* allocate(std::size_t n);

            // Deallocates a previous allocation made of size n.
            void deallocate(void* p, std::size_t n);

            // Returns if address p is from storage within the arena.
            bool is_in_arena(void* p) const noexcept;
    */
    template<std::size_t N, class T=char>
    struct arena;

    /*
        arena_alloc<T, new_fallback=true> is a C++11 conforming allocator that must be initialized
        with an arena object, and will allocate using its storage. If the arena runs out of storage
        allocations will be handled using ::operator new if new_fallback is true, else
        std::bad_alloc will be thrown.

        Example usage:

            // First 1024 bytes of memory usage by j will be served by arena:
            op::arena<1024> arena;
            std::basic_string<char, arena_alloc<int>> v{arena}; 

            // Beware reallocation and alignment issues. Use reserve immediately if you want to
            // store a particular amount of elements in a container:
            op::arena<64, double> arena;
            std::vector<int, arena_alloc<int>> v{arena}; 
            v.reserve(arena.size());

            // You can reuse the same arena for multiple allocators:
            op::arena<1024 * 8> arena; // 8 kilobytes of working memory 
            std::stack<int, arena_alloc<int>> s{arena};
            std::map<int, double, arena_alloc<std::pair<const int, double>>> m{arena};

            // And for memory constrained environments you can disallow dynamic memory allocation
            // entirely (resulting in std::bad_alloc if the arena runs out of free space):
            op::arena<1024> arena;
            std::set<double, arena_alloc<double, false>> s;
    */
    template<class T, bool new_fallback=true>
    class arena_alloc;
}



// Implementation.
namespace op {
    namespace detail {
        // Workaround, some compilers only define ::max_align_t.
        inline constexpr std::size_t min_arena_alignment() {
            using namespace std; 
            return alignof(max_align_t);
        }

        // We use a size-agnostic base class to avoid leaking the size parameter in to the
        // arena_alloc type.
        class arena_base {
        public:
            arena_base(char* storage_begin, char* storage_end) noexcept
            : storage_begin(storage_begin), storage_end(storage_end)
            , cursor(storage_begin), freelist(uint32_t(-1)) { }

            ~arena_base() noexcept { cursor = nullptr; }
            arena_base(const arena_base&) = delete;
            arena_base& operator=(const arena_base&) = delete;

            std::size_t used() const noexcept { return cursor - storage_begin; }
            void reset() noexcept { cursor = storage_begin; }

            void* allocate(std::size_t n) {
                assert(cursor != nullptr && "arena_alloc outlived arena");
                n = align(n);

                // Find first fitting block in freelist.
                uint32_t* previous_next = &freelist;
                uint32_t cur_offset = freelist;
                while (cur_offset + 1 != 0) {
                    Block* cur = get_block(cur_offset);

                    if (cur->size >= n) {
                        if (cur->size == n) {
                            *previous_next = cur->next_block;
                        } else {
                            Block *rest_block = get_block(cur_offset + n);
                            *previous_next = cur_offset + n;
                            rest_block->size = cur->size - n;
                            rest_block->next_block = cur->next_block;
                        }

                        return storage_begin + cur_offset;
                    }

                    previous_next = &cur->next_block;
                    cur_offset = cur->next_block;
                }

                // Nothing suitable found in freelist, get from storage.
                if (cursor + n <= storage_end) {
                    void* result = cursor;
                    cursor += n;
                    return result;
                }

                // No storage left.
                return nullptr;
            }

            void deallocate(void* p_, std::size_t n) {
                assert(cursor != nullptr && "arena_alloc outlived arena");
                assert(is_in_arena(p_) && "arena deallocate called with pointer outside arena");

                char* p = static_cast<char*>(p_);
                n = align(n);

                // Attempt to merge block back into storage.
                if (p + n == cursor) {
                    cursor = p;
                } else {
                    // Put block into freelist.
                    uint32_t block_offset = p - storage_begin;
                    Block* block = get_block(block_offset);

                    // No block in freelist yet or we should be first element.
                    if (block_offset < freelist) {
                        block->size = n;
                        block->next_block = freelist;
                        freelist = block_offset;
                    } else {
                        // Find biggest element smaller than us.
                        uint32_t cur_offset = freelist;
                        Block* cur = get_block(cur_offset);
                        while (block_offset > cur->next_block) {
                            cur_offset = cur->next_block;
                            cur = get_block(cur_offset);
                        }

                        // Merge into preceding block or insert.
                        if (cur_offset + cur->size == block_offset) {
                            cur->size += n;
                            block_offset = cur_offset;
                            block = cur;
                        } else {
                            block->size = n;
                            block->next_block = cur->next_block;
                            cur->next_block = block_offset;
                        }
                    }

                    // Merge this block with adjacent block to the right if possible.
                    if (block_offset + block->size == block->next_block) {
                        Block* next_block = get_block(block->next_block);
                        block->size += next_block->size;
                        block->next_block = next_block->next_block;
                    }
                }
            }

            bool is_in_arena(void* p_) const noexcept {
                char* p = static_cast<char*>(p_);
                return storage_begin <= p && p < storage_end;
            }

        private:
            struct Block {
                uint32_t next_block;
                uint32_t size;
            };

            Block* get_block(uint32_t block_offset) {
                return reinterpret_cast<Block*>(storage_begin + block_offset);
            }

            char* storage_begin;
            char* storage_end;
            char* cursor;
            uint32_t freelist;

        protected:
            static constexpr std::size_t alignment =
                min_arena_alignment() < sizeof(Block) ? sizeof(Block) : min_arena_alignment();
            
            static constexpr std::size_t align(std::size_t n) {
                return ((n + alignment - 1) / alignment) * alignment;
            }
        };

        template<class T, bool new_fallback>
        class arena_alloc_fallback_base {
        protected:
            T* fallback_allocate(std::size_t n) {
                return static_cast<T*>(::operator new(n*sizeof(T)));
            }

            void fallback_deallocate(T* p, std::size_t n) noexcept {
                ::operator delete(p);
            }
        };

        template<class T>
        class arena_alloc_fallback_base<T, false> {
        protected:
            T* fallback_allocate(std::size_t n) {
                throw std::bad_alloc();
            }

            void fallback_deallocate(T* p, std::size_t n) noexcept {
                assert(false && "arena_alloc deallocate called with pointer not from arena");
            }
        };
    }

    template<std::size_t N, class T>
    struct arena : detail::arena_base {
        static_assert(N * sizeof(T) < (1ull << 32) - 8, "arena can only handle up to 2^32-8 bytes");
        arena() noexcept : detail::arena_base(storage, storage + N) { }
        static constexpr std::size_t size() noexcept { return N; }
        alignas(alignment) char storage[align(N * sizeof(T))];
    };


    template<class T, bool new_fallback>
    class arena_alloc : private detail::arena_alloc_fallback_base<T, new_fallback> {
    public:
        typedef T value_type;
        template<class U> struct rebind { typedef arena_alloc<U, new_fallback> other; };

        arena_alloc(detail::arena_base& a) : a(a) { }
        template<class U, bool Un>
        arena_alloc(const arena_alloc<U, Un>& alloc) noexcept : a(alloc.a) { }
        arena_alloc(const arena_alloc&) = default;
        arena_alloc& operator=(const arena_alloc&) = delete;

        T* allocate(std::size_t n) {
            T* result = static_cast<T*>(a.allocate(n*sizeof(T)));
            if (result) return result;
            return this->fallback_allocate(n);
        }

        void deallocate(T* p, std::size_t n) noexcept {
            if (a.is_in_arena(p)) a.deallocate(p, n*sizeof(T));
            else this->fallback_deallocate(p, n);
        }

        template<class T_, bool Tn, class U, bool Un>
        friend bool operator==(const arena_alloc<T_, Tn>& x, const arena_alloc<U, Un>& y) noexcept;

    private:
        detail::arena_base& a;
    };

    template<class T, bool Tn, class U, bool Un>
    inline bool operator==(const arena_alloc<T, Tn>& x, const arena_alloc<U, Un>& y) noexcept {
        return &x.a == &y.a;
    }

    template<class T, bool Tn, class U, bool Un>
    inline bool operator!=(const arena_alloc<T, Tn>& x, const arena_alloc<U, Un>& y) noexcept {
        return !(x == y);
    }
}

#endif
