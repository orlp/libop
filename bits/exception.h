#ifndef OP_EXCEPTION_H
#define OP_EXCEPTION_H

#include <stdexcept>
#include <memory>
#include <string>

namespace op {
    class BaseException : public virtual std::exception {
    public:
        explicit BaseException(const std::string& msg) : msg_storage(std::make_shared<std::string>(msg)) { }
        explicit BaseException(const char* msg) : msg_storage(std::make_shared<std::string>(msg)) { }
        
        virtual const char* what() const noexcept { return msg_storage->c_str(); }
        
    private:
        // shared_ptr to make copy constructor noexcept
        std::shared_ptr<std::string> msg_storage;
    };
}

#endif
