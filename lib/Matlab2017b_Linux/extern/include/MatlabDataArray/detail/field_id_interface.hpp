/* Copyright 2016 The MathWorks, Inc. */

#ifndef FIELD_ID_INTERFACE_HPP_
#define FIELD_ID_INTERFACE_HPP_

#include "publish_util.hpp"

namespace matlab {
    namespace data {
        namespace detail {
            class RefCounted;
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::RefCounted* field_id_create_impl(const char* str, size_t len);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void field_id_destroy_impl(matlab::data::detail::RefCounted* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool field_id_is_equal(matlab::data::detail::RefCounted* lhs, matlab::data::detail::RefCounted* rhs);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void field_id_get_string(matlab::data::detail::RefCounted* impl, const char** str, size_t* len);


#endif
