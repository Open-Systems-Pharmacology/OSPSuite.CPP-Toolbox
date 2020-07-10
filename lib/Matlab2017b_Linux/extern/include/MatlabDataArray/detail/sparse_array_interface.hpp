/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef SPARSE_ARRAY_INTERFACE_HPP_
#define SPARSE_ARRAY_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace impl {
            class ArrayImpl;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void sparse_array_get_num_nonzero_elements(matlab::data::impl::ArrayImpl* impl, size_t * val);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void sparse_array_get_index(matlab::data::impl::ArrayImpl* impl, matlab::data::detail::IteratorImpl* itImpl, size_t * row, size_t * col);

#endif
