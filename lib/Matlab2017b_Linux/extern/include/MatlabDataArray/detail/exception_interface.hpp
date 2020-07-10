/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef EXCEPTION_INTERFACE_HPP_
#define EXCEPTION_INTERFACE_HPP_

#include "publish_util.hpp"
#include <exception>
#include <string>

namespace matlab {
    namespace data {
        namespace detail {
            template <typename Base, matlab::data::ExceptionType Type>
            class ArrayException : public Base, public std::exception {
              public:
                ArrayException(std::string str) MW_NOEXCEPT :
                    fStringTxt(std::move(str)) {}

                ~ArrayException() MW_NOEXCEPT {}

                ArrayException(const ArrayException& other) MW_NOEXCEPT :
                    fStringTxt(other.fStringTxt) {}

                const char* what() const MW_NOEXCEPT {
                    return fStringTxt.c_str();
                }
                
              private:
                std::string fStringTxt;
            };
        }
    }
}
                
#endif
