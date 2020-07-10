/* Copyright 2015-2017 The MathWorks, Inc. */

#ifndef REFERENCE_INTERFACE_HPP_
#define REFERENCE_INTERFACE_HPP_

#include "publish_util.hpp"
#include <stdint.h>
#include <cstddef>

namespace matlab {
    namespace data {
        class Array;

        namespace detail {
            class ReferenceImpl;
            class IteratorImpl;
            class ForwardIteratorImpl;
        }
        namespace impl {
            class ArrayImpl;
        }
    }
}

// Start: Untyped APIs - need to validate type
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_add_index(matlab::data::detail::ReferenceImpl* impl,
                                                              size_t idx);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_add_string_index(matlab::data::detail::ReferenceImpl* impl,
                                                                     const char* idx,
                                                                     size_t len);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_reference_value(matlab::data::detail::ReferenceImpl* impl,
                                                                        matlab::data::impl::ArrayImpl* rhs);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void reference_get_reference_value(matlab::data::detail::ReferenceImpl* impl,
                                                                bool unshareData,
                                                                matlab::data::detail::ReferenceImpl** retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void reference_destroy_impl(matlab::data::detail::ReferenceImpl* impl);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_pod_value(matlab::data::detail::ReferenceImpl* impl,
                                                                   void* rhs);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_complex_value(matlab::data::detail::ReferenceImpl* impl,
                                                              void* real,
                                                              void* imag );

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_get_pod_value(matlab::data::detail::ReferenceImpl* impl,
                                                                   void** retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_get_complex_value(matlab::data::detail::ReferenceImpl* impl, 
                                                                       void** real,  
                                                                       void** imag ); 

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_validate_index(matlab::data::detail::ReferenceImpl* impl);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_validate_type(matlab::data::detail::ReferenceImpl* impl,
                                                                  int type);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_array_reference_validate_type(matlab::data::detail::ReferenceImpl* impl,
                                                                              int type,
                                                                              bool unshare);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_get_container_type(matlab::data::detail::ReferenceImpl* impl, int *type);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_string(matlab::data::detail::ReferenceImpl* impl, const char* val, size_t len);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_char16_string(matlab::data::detail::ReferenceImpl* impl, const char16_t* val, size_t len);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int reference_set_missing_char16_string(matlab::data::detail::ReferenceImpl* impl);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int char_array_ref_get_string(matlab::data::detail::ReferenceImpl* impl, char16_t const** str, size_t* strLen);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int char_array_ref_get_ascii(matlab::data::detail::ReferenceImpl* impl, char16_t const** str, size_t* strLen);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void enum_ref_get_class(matlab::data::detail::ReferenceImpl* impl,
                                                              char** str,
                                                              size_t* strlen);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::ForwardIteratorImpl* struct_array_ref_begin_id(matlab::data::detail::ReferenceImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::ForwardIteratorImpl* struct_array_ref_end_id(matlab::data::detail::ReferenceImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C size_t struct_array_ref_get_num_fields(matlab::data::detail::ReferenceImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void sparse_array_ref_get_num_nonzero_elements(matlab::data::detail::ReferenceImpl* impl, size_t * val);

// End Untyped APIs


// Start: Array specific APIs:

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_reference_get_dimensions(matlab::data::detail::ReferenceImpl* impl,
                                                                          size_t* numDims,
                                                                          size_t** dims);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_reference_get_type(matlab::data::detail::ReferenceImpl* impl, 
                                                                    int* retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_reference_get_num_elements(matlab::data::detail::ReferenceImpl* impl, size_t* retVal);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_reference_is_empty(matlab::data::detail::ReferenceImpl* impl, bool* retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::IteratorImpl* array_reference_begin(matlab::data::detail::ReferenceImpl* impl, bool unshare);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::IteratorImpl* array_reference_end(matlab::data::detail::ReferenceImpl* impl, bool unshare);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int array_reference_get_reference(matlab::data::detail::ReferenceImpl* impl,
                                                                        bool unshareData,
                                                                        size_t idx,
                                                                        matlab::data::detail::ReferenceImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void array_reference_set_value(matlab::data::detail::ReferenceImpl* impl,
                                                                     matlab::data::impl::ArrayImpl* rhs);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void array_reference_shared_copy(matlab::data::detail::ReferenceImpl* impl,
                                                                       matlab::data::impl::ArrayImpl** retVal);


// End: Array specific APIs:


// Start: Typed APIs:

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_reference_set_pod_value(matlab::data::detail::ReferenceImpl* impl,
                                                                        int type,
                                                                        void* rhs);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_reference_set_complex_value(matlab::data::detail::ReferenceImpl* impl,
                                                                            int type,
                                                                            void* real,
                                                                            void* imag );

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_reference_get_pod_value(matlab::data::detail::ReferenceImpl* impl,
                                                                         void** retVal);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int typed_reference_get_complex_value(matlab::data::detail::ReferenceImpl* impl, 
                                                                             void** real,  
                                                                             void** imag ); 
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool reference_matlab_string_equal_string(matlab::data::detail::ReferenceImpl* impl, const char* str, size_t strLen);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool reference_matlab_string_equal_string16(matlab::data::detail::ReferenceImpl* impl, const char16_t* str, size_t strLen);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C bool reference_matlab_string_equal_reference(matlab::data::detail::ReferenceImpl* lhs, matlab::data::detail::ReferenceImpl* rhs);

// End: Typed APIs:

#endif
