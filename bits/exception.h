#ifndef OP_EXCEPTION_H
#define OP_EXCEPTION_H

#include <algorithm>
#include <exception>
#include <memory>
#include <string>

/*
   A base exception that stores a message and plays nicely with virtual inheritance.

    In order to make easy and proper exception hierarchies with this base class, follow the
    following rules:
    
     1. Every exception inherits publicly and virtually from its parent exceptions.
    
     2. Every exception declares a protected default constructor and defines a protected
        constructor initializing all data members.
    
     3. Every exception that is supposed to be Constructible defines a public constructor that
        directly calls the constructors defined in 2 for every ancestor.
    
     4. All copy constructors should be noexcept.
*/

namespace op {
    struct BaseException : virtual std::exception {
        explicit BaseException(std::string msg)
        : msg_storage(std::make_shared<std::string>(std::move(msg))) { }
        
        virtual const char* what() const noexcept { return msg_storage->c_str(); }

    protected:
        BaseException();

    private:
        // shared_ptr to make copy constructor noexcept.
        std::shared_ptr<std::string> msg_storage;
    };
}

#endif
