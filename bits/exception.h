#ifndef OP_EXCEPTION_H
#define OP_EXCEPTION_H

#include <exception>
#include <memory>
#include <string>


namespace op {
    // A base exception that stores a message and plays nicely with virtual inheritance.
    class BaseException : public virtual std::exception {
    public:
        explicit BaseException(std::string msg)
        : msg_storage(std::make_shared<std::string>(std::move(msg))) { }
        
        virtual const char* what() const noexcept { return msg_storage->c_str(); }
        
    private:
        // Shared_ptr to make copy constructor noexcept.
        std::shared_ptr<std::string> msg_storage;
    };
}

#endif
