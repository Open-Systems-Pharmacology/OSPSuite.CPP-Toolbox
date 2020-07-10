/* Copyright 2015 The MathWorks, Inc. */

#ifndef STRING_INTERFACE_HPP_
#define STRING_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>
#include <string>
#include <algorithm>

namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;

            inline bool isAscii7(const String &str) {
                return (std::find_if(str.begin(), str.end(), [](CHAR16_T ch) {return (ch & 0xFF80) != 0; }) == str.end());
            }

            inline bool isAscii7(const std::string &str) {
                return (std::find_if(str.begin(), str.end(), [](char ch) {return (ch & 0x80) != 0; }) == str.end());
            }
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int string_get_value(matlab::data::detail::ReferenceImpl* impl,
    char16_t** str,
    size_t* strlen);

#endif
