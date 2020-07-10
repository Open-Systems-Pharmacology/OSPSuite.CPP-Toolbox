/* Copyright 2017 The MathWorks, Inc. */

#ifndef STRUCT_INTERFACE_HPP_
#define STRUCT_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;
            class ForwardIteratorImpl;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int struct_get_field_copy(matlab::data::detail::ReferenceImpl* impl,
                                                                const char* stringIndex,
                                                                size_t stringIndexLen,
                                                                matlab::data::impl::ArrayImpl** retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int struct_reference_get_index(matlab::data::detail::ReferenceImpl* impl,
                                                                     const char* stringIndex,
                                                                     size_t stringIndexLen,
                                                                     size_t* retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::ForwardIteratorImpl* struct_array_begin_id(matlab::data::impl::ArrayImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::ForwardIteratorImpl* struct_array_end_id(matlab::data::impl::ArrayImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C size_t struct_array_get_num_fields(matlab::data::impl::ArrayImpl* impl);


#endif
