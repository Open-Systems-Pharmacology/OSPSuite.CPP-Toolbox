/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef ARRAY_FACTORY_INTERFACE_HPP_
#define ARRAY_FACTORY_INTERFACE_HPP_

#include "publish_util.hpp"
#include "ExceptionHelpers.hpp"

namespace matlab {
    namespace data {
        namespace impl {
            class ArrayFactoryImpl;
            class ArrayImpl;
        }
        namespace detail {
            class NameListImpl;

            
        }
    }
}

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void array_factory_destroy_impl(matlab::data::impl::ArrayFactoryImpl* impl);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::impl::ArrayFactoryImpl* array_factory_create();

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_array_with_dims(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                 int arrayType,
                                                                 size_t* dims, 
                                                                 size_t numDims,
                                                                 matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_array_with_dims_and_data(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                          int arrayType,
                                                                          size_t* dims, 
                                                                          size_t numDims,
                                                                          const void* const dataStart,
                                                                          size_t numEl,
                                                                          matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_scalar_array(matlab::data::impl::ArrayFactoryImpl* impl,
                                                              int arrayType,
                                                              const void* data,
                                                              matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_scalar_string(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                               const char16_t* data,
                                                               size_t strlen,
                                                               matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_scalar_missing_string(matlab::data::impl::ArrayFactoryImpl* impl,
                                                                       matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_char_array_from_char16_string(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                               const char16_t* data,
                                                                               size_t strlen,
                                                                               matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C  int create_char_array_from_string(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                     const char* data,
                                                                     size_t strlen,
                                                                     matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_buffer(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                        void** buffer, void(**deleter)(void*), int dataType, size_t numElements);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_array_from_buffer(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                   int arrayType,
                                                                   size_t* dims, 
                                                                   size_t numDims,
                                                                   void* buffer,
                                                                   void(* deleter)(void*),
                                                                   matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_sparse_array_from_buffer(matlab::data::impl::ArrayFactoryImpl* impl, 
                                                                          int arrayType,
                                                                          size_t* dims, 
                                                                          size_t numDims,
                                                                          size_t nnz,
                                                                          void* dataBuffer,
                                                                          void(* dataDeleter)(void*),
                                                                          size_t* rowsBuffer,
                                                                          void(* rowsDeleter)(void*),
                                                                          size_t* colsBuffer,
                                                                          void(* colsDeleter)(void*),
                                                                          matlab::data::impl::ArrayImpl**);

MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C matlab::data::detail::NameListImpl* create_names(size_t num);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void add_name(matlab::data::detail::NameListImpl* impl, const char* name, size_t nameLen);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_struct_array(matlab::data::impl::ArrayFactoryImpl* impl, size_t* dims, size_t numDims, matlab::data::detail::NameListImpl* names, matlab::data::impl::ArrayImpl**);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C void names_destroy_impl(matlab::data::detail::NameListImpl* impl);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_enum_array_with_enums(matlab::data::impl::ArrayFactoryImpl* impl, size_t* dims, size_t numDims, const char* cls, size_t clsSize, matlab::data::detail::NameListImpl* names, matlab::data::impl::ArrayImpl** retVal);
MATLAB_EXTDATA_PUBLISHED_API_EXTERN_C int create_enum_array(matlab::data::impl::ArrayFactoryImpl* impl, size_t* dims, size_t numDims, const char* cls, size_t clsSize, matlab::data::impl::ArrayImpl** retVal);


namespace matlab {
    namespace data {
        namespace detail {
            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
            typename std::enable_if<std::is_arithmetic<T>::value && !std::is_same<T, bool>::value, TypedArray<T>>::type 
            createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {
                matlab::data::impl::ArrayImpl* aImpl = nullptr;
                throwIfError(create_array_with_dims_and_data(
                                 impl,
                                 static_cast<int>(GetArrayType<T>::type),
                                 &dims[0], 
                                 dims.size(),
                                 &(*(begin)),
                                 (end - begin),
                                 &aImpl));
                return matlab::data::detail::Access::createObj<TypedArray<T>>(aImpl);
            }

            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
            typename std::enable_if<std::is_same<T, bool>::value, TypedArray<bool>>::type  
            createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {
                matlab::data::impl::ArrayImpl* aImpl = nullptr;
                throwIfError(create_array_with_dims(
                                 impl,
                                 static_cast<int>(GetArrayType<T>::type),
                                 &dims[0],
                                 dims.size(),
                                 &aImpl));
                auto retVal = matlab::data::detail::Access::createObj<TypedArray<T>>(aImpl);

                auto it = begin;
                for (auto& elem : retVal) {
                    elem = *it;
                    if (++it == end) {
                        break;
                    }
                }
                return retVal;
            }

            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
            typename std::enable_if<matlab::data::is_complex<T>::value || std::is_same<T, String>::value, TypedArray<typename GetReturnType<T>::type>>::type  
            createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {
                matlab::data::impl::ArrayImpl* aImpl = nullptr;
                throwIfError(create_array_with_dims(
                                 impl,
                                 static_cast<int>(GetArrayType<typename GetReturnType<T>::type>::type),
                                 &dims[0],
                                 dims.size(),
                                 &aImpl));
                auto retVal = matlab::data::detail::Access::createObj<TypedArray<typename GetReturnType<T>::type>>(aImpl);

                auto it = begin;
                for (auto elem : retVal) {
                    elem = *it;
                    if (++it == end) {
                        break;
                    }
                }
                return retVal;
            }

            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
            typename std::enable_if<std::is_same<T, Array>::value, TypedArray<T>>::type  
            createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {
                matlab::data::impl::ArrayImpl* aImpl = nullptr;
                throwIfError(create_array_with_dims(
                                 impl,
                                 static_cast<int>(GetArrayType<T>::type),
                                 &dims[0],
                                 dims.size(),
                                 &aImpl));
                             
                auto retVal = matlab::data::detail::Access::createObj<TypedArray<T>>(aImpl);

                auto it = begin;
                for (auto elem : retVal) {
                    elem = T(*it);
                    if (++it == end) {
                        break;
                    }
                }
                return retVal;
            }


            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
                      typename std::enable_if<std::is_same<T, const char *>::value || std::is_same<T, std::string>::value, TypedArray<MATLABString>>::type
            createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {                
                if (std::find_if(begin, end, [](std::string str) {return (detail::isAscii7(str) == false); }) != end) {
                    throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                }

		matlab::data::impl::ArrayImpl* aImpl = nullptr;
		throwIfError(create_array_with_dims(
		             impl,
			     static_cast<int>(ArrayType::MATLAB_STRING),
			     &dims[0],
			     dims.size(),
			     &aImpl));
		auto retVal = matlab::data::detail::Access::createObj<TypedArray<MATLABString>>(aImpl);

		auto it = begin;
		for (auto elem : retVal) {
		    elem = *it;
		    if (++it == end) {
     		        break;
		    }
		}
		return retVal;
	    }

            template <typename ItType,
                typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
                typename std::enable_if<std::is_same<T, MATLABString>::value, TypedArray<MATLABString>>::type
                createArrayWithIterator(matlab::data::impl::ArrayFactoryImpl* impl, ArrayDimensions dims, ItType begin, ItType end) {

                matlab::data::impl::ArrayImpl* aImpl = nullptr;
                throwIfError(create_array_with_dims(
                    impl,
                    static_cast<int>(ArrayType::MATLAB_STRING),
                    &dims[0],
                    dims.size(),
                    &aImpl));
                auto retVal = matlab::data::detail::Access::createObj<TypedArray<MATLABString>>(aImpl);

                auto it = begin;
                for (auto elem : retVal) {
                    MATLABString x = *it;
                    if (x) {
                        elem = x;
                    }
                    if (++it == end) {
                        break;
                    }
                }
                return retVal;
            }
        }
    }
}


#endif
