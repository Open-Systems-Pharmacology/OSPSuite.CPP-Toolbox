/* Copyright 2015-2017 The MathWorks, Inc. */

#ifndef TYPED_ITERATOR_INTERFACE_HPP_
#define TYPED_ITERATOR_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace detail {
            class IteratorImpl;
            class ReferenceImpl;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_destroy_impl(matlab::data::detail::IteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::IteratorImpl* typed_iterator_clone(matlab::data::detail::IteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_plus_plus(matlab::data::detail::IteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_minus_minus(matlab::data::detail::IteratorImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_increment(matlab::data::detail::IteratorImpl* impl, ptrdiff_t);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_decrement(matlab::data::detail::IteratorImpl* impl, ptrdiff_t);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool typed_iterator_equal(matlab::data::detail::IteratorImpl* impl, matlab::data::detail::IteratorImpl* rhs);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C ptrdiff_t typed_iterator_distance_to(matlab::data::detail::IteratorImpl* impl, matlab::data::detail::IteratorImpl* rhs);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_get_pod_value(matlab::data::detail::IteratorImpl* impl, void** val);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void typed_iterator_get_proxy(matlab::data::detail::IteratorImpl* impl, matlab::data::detail::ReferenceImpl** val);

#endif
