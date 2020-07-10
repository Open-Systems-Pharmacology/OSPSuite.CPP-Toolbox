/* Copyright 2015 The MathWorks, Inc. */

#ifndef FORWARD_ITERATOR_INTERFACE_HPP_
#define FORWARD_ITERATOR_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace detail {
            class ForwardIteratorImpl;
            class RefCounted;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void forward_iterator_destroy_impl(matlab::data::detail::ForwardIteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::ForwardIteratorImpl* forward_iterator_clone(matlab::data::detail::ForwardIteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void forward_iterator_plus_plus(matlab::data::detail::ForwardIteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool forward_iterator_equal(matlab::data::detail::ForwardIteratorImpl* impl, matlab::data::detail::ForwardIteratorImpl* rhs);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void forward_iterator_get_ref(matlab::data::detail::ForwardIteratorImpl* impl, matlab::data::detail::RefCounted** val);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void forward_iterator_increment(matlab::data::detail::ForwardIteratorImpl* impl, size_t d);
#endif
