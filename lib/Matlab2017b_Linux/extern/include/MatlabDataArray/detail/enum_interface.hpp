/* Copyright 2015 The MathWorks, Inc. */

#ifndef ENUM_INTERFACE_HPP_
#define ENUM_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;
        }
        namespace impl {
            class ArrayImpl;
        }
    }
}


MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void enum_get_class(matlab::data::impl::ArrayImpl* impl,
                                                          char** str,
                                                          size_t* strlen);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int enum_get_value(matlab::data::detail::ReferenceImpl* impl,
                                                          char** str,
                                                          size_t* strlen);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int enum_set_value(matlab::data::detail::ReferenceImpl* impl,
                                                         const char* str,
                                                         size_t strlen);

#endif
