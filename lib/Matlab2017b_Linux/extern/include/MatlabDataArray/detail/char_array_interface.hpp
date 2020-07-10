/* Copyright 2016 The MathWorks, Inc. */

#ifndef CHAR_ARRAY_INTERFACE_HPP_
#define CHAR_ARRAY_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace impl {
            class ArrayImpl;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void char_array_get_string(matlab::data::impl::ArrayImpl* impl, char16_t const** str, size_t* strLen);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int char_array_get_ascii(matlab::data::impl::ArrayImpl* impl, char16_t const** str, size_t* strLen);

#endif
