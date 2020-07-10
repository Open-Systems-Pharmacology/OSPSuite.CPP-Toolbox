/* Copyright 2015-2017 The MathWorks, Inc. */

#ifndef ARRAY_INTERFACE_HPP_
#define ARRAY_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>

namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;
            class IteratorImpl;
        }
        namespace impl {
            class ArrayImpl;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void array_destroy_impl(matlab::data::impl::ArrayImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::impl::ArrayImpl* array_create_empty();
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_get_type(matlab::data::impl::ArrayImpl* impl, int *type); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool array_is_empty(matlab::data::impl::ArrayImpl* impl); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C size_t array_get_number_of_elements(matlab::data::impl::ArrayImpl* impl); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void array_get_dimensions(matlab::data::impl::ArrayImpl* impl,
                                                                size_t* numDims,
                                                                size_t** dims);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool array_unshare(matlab::data::impl::ArrayImpl* impl, bool isUserArrayUnique, matlab::data::impl::ArrayImpl** newImpl); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool array_is_empty(matlab::data::impl::ArrayImpl* impl); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_create_reference(matlab::data::impl::ArrayImpl* impl,
                                                                 size_t idx,
                                                                 matlab::data::detail::ReferenceImpl** retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::IteratorImpl* typed_array_begin(matlab::data::impl::ArrayImpl* impl, bool unshare);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::IteratorImpl* typed_array_end(matlab::data::impl::ArrayImpl* impl, bool unshare);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_array_is_valid_conversion(int lhsDataType, int rhsDataType, bool *result);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_release_buffer(matlab::data::impl::ArrayImpl* impl, bool isUserArrayUnique, void** buffer, void(**deleter)(void*));

#endif
