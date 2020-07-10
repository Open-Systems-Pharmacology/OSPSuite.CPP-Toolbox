/* Copyright 2014-2017 The MathWorks, Inc. */

#ifndef MATLAB_DATA_ARRAY_FACTORY_HPP_
#define MATLAB_DATA_ARRAY_FACTORY_HPP_

#include "TypedArray.hpp"
#include "ArrayDimensions.hpp"
#include "CharArray.hpp"
#include "StructArray.hpp"
#include "SparseArray.hpp"
#include "Enumeration.hpp"
#include "EnumArray.hpp"
#include "Exception.hpp"
#include "GetArrayType.hpp"
#include "String.hpp"
#include "matlab_extdata_defs.hpp"
#include "TypedIterator.hpp"
#include "Optional.hpp"
#include "GetReturnType.hpp"

#include "detail/array_factory_interface.hpp"
#include "detail/string_interface.hpp"
#include "detail/HelperFunctions.hpp"
#include "detail/ExceptionHelpers.hpp"
#include "detail/publish_util.hpp"

#include <string>
#include <initializer_list>
#include <stdint.h>
#include <iterator>
#include <vector>

namespace matlab {
    namespace data {
        namespace impl {
            class ArrayFactoryImpl;
        }

        namespace detail {
            class NameListImpl;
        }
        

        /**
         * The ArrayFactory base class provides all of the general APIs that 
         * each of the concrete factories need to support.
         */
        class ArrayFactory {
          public:


            /**
             * Construct an ArrayFactory
             *
             * @throw FailedToLoadLibMatlabDataArrayException if necessary libraries are not found
             */
            ArrayFactory() :
                pImpl(array_factory_create(), 
                      [](impl::ArrayFactoryImpl* ptr) {
                          array_factory_destroy_impl(ptr);
                      }) {
                detail::throwIfError(pImpl.get() ? 0 : static_cast<int>(ExceptionType::FailedToLoadLibMatlabDataArray));
            }
                            
            /**
             * ArrayFactory destructor
             *
             * @throw none
             */
            ~ArrayFactory() MW_NOEXCEPT {}

            /**
             * Creates a TypedArray<T> with the given dimensions
             *
             * @param dims - the dimensions for the Array
             *
             * @return TypedArray<T> - an Array with the appropriate type and dimensions
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw InvalidArrayType - if array type of the new array is not recognized as valid
             */
            template<typename T> 
            TypedArray<T> createArray(ArrayDimensions dims) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_array_with_dims(
                                         pImpl.get(),
                                         static_cast<int>(GetArrayType<T>::type),
                                         &dims[0],
                                         dims.size(),
                                         &impl));
                return matlab::data::detail::Access::createObj<TypedArray<T>>(impl);
            }

            /**
             * Creates a TypedArray<T> with the given dimensions and filled in with the
             * data specified by the begin and end iterators. Data is copied and must be in column major order.
             * The data type for the
             * array is determined by the value_type of the iterator.
             *
             * @param dims - the dimensions for the Array
             * @param begin - start of the user supplied data
             * @param end - end of the user supplied data
             *
             * @return TypedArray<T> - an Array with the appropriate type, dimensions and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            template <typename ItType,
                      typename T = typename std::remove_cv<typename std::iterator_traits<ItType>::value_type>::type>
            TypedArray< typename GetReturnType<T>::type> createArray(ArrayDimensions dims, ItType begin,  ItType end) {
                return detail::createArrayWithIterator(pImpl.get(), std::move(dims), begin, end);
            }

            /**
             * Creates a TypedArray<T> with the given dimensions and filled in with the
             * supplied data in the initializer list. Data must be in column major order.
             *
             * @param dims - the dimensions for the Array
             * @param data - initializer list containing the data
             *
             * @return TypedArray<T> - an Array with the appropriate type, dimensions and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */

            template<typename T> 
            TypedArray<typename GetReturnType<T>::type> createArray(ArrayDimensions dims,
                                      std::initializer_list<T> data) {
                return detail::createArrayWithIterator(pImpl.get(), std::move(dims), data.begin(), data.end());
            };

            /**
             * Creates a TypedArray<T> with the given dimensions and filled in with the
             * data specified by C-style begin and end pointers. Data is copied and must be in column major order.
             * This methods supports all arithmetic types. 
             *
             * @param dims - the dimensions for the Array
             * @param T*  - start of the user supplied data
             * @param T*  - end of the user supplied data
             *
             * @return TypedArray<T> - an Array with the appropriate type, dimensions and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw InvalidArrayType - if array type of the new array is not recognized as valid
             */
            template <typename T>
            typename std::enable_if<std::is_arithmetic<T>::value, TypedArray<T>>::type
            createArray(ArrayDimensions dims, const T* const begin, const T* const end) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_array_with_dims_and_data(
                                         pImpl.get(),
                                         static_cast<int>(GetArrayType<T>::type),
                                         &dims[0],
                                         dims.size(),
                                         begin,
                                         (end - begin),
                                         &impl));
                return matlab::data::detail::Access::createObj<TypedArray<T>>(impl);
            }

            /**
             * Creates a TypedArray<T> with the given dimensions and filled in with the
             * data specified by C-style begin and end pointers. Data is copied and must be in column major order.
             * This methods supports std::complex types.
             *
             * @param dims - the dimensions for the Array
             * @param T*  - start of the user supplied data
             * @param T*  - end of the user supplied data
             *
             * @return TypedArray<T> - an Array with the appropriate type, dimensions and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            template <typename T>
            typename std::enable_if<matlab::data::is_complex<T>::value, TypedArray<T>>::type
            createArray(ArrayDimensions dims, const T* const begin, const T* const end) {
                auto retVal = createArray<T>(std::move(dims));

                auto it = begin;
                for (auto elem : retVal) {
                    elem = *it;
                    if (++it == end) {
                        break;
                    }
                }
                return retVal;
            }

            /**
             * Creates a scalar TypedArray<T> with the given value. This methods supports
             * all arithmetic types.
             *
             * @param val - the value to be inserted into the scalar
             *
             * @return TypedArray<T> - an Array with the appropriate type and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw InvalidArrayTypeException - if array type is not recognized as valid
             */
            template<typename T>
            typename std::enable_if<std::is_arithmetic<T>::value, TypedArray<T>>::type 
            createScalar(const T val) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_scalar_array(pImpl.get(),
                                                         static_cast<int>(GetArrayType<T>::type), 
                                                         &val,
                                                         &impl));
                return matlab::data::detail::Access::createObj<TypedArray<T>>(impl);
            }

            /**
             * Creates a scalar TypedArray<T> with the given value. This methods supports
             * all complex types.
             *
             * @param val - the value to be inserted into the scalar
             *
             * @return TypedArray<T> - an Array with the appropriate type and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            template<typename T> 
            typename std::enable_if<matlab::data::is_complex<T>::value, TypedArray<T>>::type createScalar(const T val) {
                auto arr = createArray<T>({1,1});
                arr[0] = val;
                return arr;
            }

            /**
             * Creates a scalar TypedArray<MATLABString> with the given value
             *
             * @param val - the string to be inserted into the scalar
             *
             * @return TypedArray<MATLABString> - a newly created TypedArray of MATLABString.
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            TypedArray<MATLABString> createScalar(const String val) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_scalar_string(pImpl.get(),
                                                          val.c_str(),
                                                          val.size(),
                                                          &impl));
                return matlab::data::detail::Access::createObj<TypedArray<MATLABString>>(impl);
            }

            /**
             * Creates a scalar TypedArray<MATLABString> with the given value
             *
             * @param val - the value to be inserted into the scalar
             *
             * @return TypedArray<MATLABString> - an Array containing data of MATLABString type
             *
             * @throw OutOfMemoryException - if the array could not be allocated
             */
            TypedArray<MATLABString> createScalar(const MATLABString val) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                if (val) {
                    String str = *val;
                    detail::throwIfError(create_scalar_string(pImpl.get(),
                                                              str.c_str(),
                                                              str.size(),
                                                              &impl));
                } else {
                    detail::throwIfError(create_scalar_missing_string(pImpl.get(), &impl));
                }
                return matlab::data::detail::Access::createObj<TypedArray<MATLABString>>(impl);
            }

           /**
            * Creates a CellArray with specified dimensions.
            *
            * @param dims - the dimensions of the CellArray
            *
            * @return CellArray - newly created Cell Array
            *
            * @throw OutOfMemoryException - if the array could not be allocated
            */
            CellArray createCellArray(ArrayDimensions dims) {
                return createArray<Array>(std::move(dims));
            }

            /**
             * Creates a CellArray with the data specified. Data is specified in column major order
             *
             * @param dims - the dimensions of the CellArray
             * @params data - elements to be inserted into the cell array. Primitive types, strings and 
             *                Arrays are allowable inputs.
             *
             * @return TypedArray<Array> - a Cell Array
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            template <typename ...Targs>
            TypedArray<Array> createCellArray(ArrayDimensions dims, Targs... data) {
                auto retVal = createArray<Array>(std::move(dims));
                setCellValues(retVal.begin(), retVal.end(), data...);
                return retVal;
            }

            /**
             * Creates a scalar TypedArray<MATLABString> from 7-bit ascii input data
             * Input is converted to UTF16 prior to creating a TypedArray<MATLABString>
             *
             * @param val - std::string to be inserted into the scalar StringArray
             *
             * @return TypedArray<MATLABString> - an Array with the appropriate type and data
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw NonAsciiCharInInputDataException - if input contains non-ascii data
             */
            TypedArray<MATLABString> createScalar(const std::string val) {
                if (!detail::isAscii7(val)) {
                    throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                }
                return createScalar(String(val.begin(), val.end())); 
            }

            /**
             * Creates a 1xn CharArray from the specified char16_t string, where n is the 
             * string length
             *
             * @param str - String containing the data to be filled into the Array
             *
             * @return CharArray - a CharArray object containing data copied from str
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */
            CharArray createCharArray(String str) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_char_array_from_char16_string(pImpl.get(), str.c_str(), str.size(), &impl));
                return matlab::data::detail::Access::createObj<CharArray>(impl);
            }
            
            /**
             * Creates a 1xn CharArray from the specified std::string, where n is the 
             * string length.  Data is checked for 7-bit ascii and converted from 
             * std::string<char> to std::string<CHAR16_T>.
             *
             * @param str - std::string containing the data to be filled into the Array
             *
             * @return CharArray - a CharArray object containing data copied from str
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw NonAsciiCharInInputDataException - if data contains non-ascii characters
             */
            CharArray createCharArray(std::string str) {
                if (!detail::isAscii7(str)) {
                    throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                }
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_char_array_from_string(pImpl.get(), str.c_str(), str.size(), &impl));
                return matlab::data::detail::Access::createObj<CharArray>(impl);
            }

            /**
             * Creates a StructArray with the given dimensions and fieldnames
             *
             * @param dims - the dimensions for the Array
             * @param fieldNames - vector of the field names for the struct
             *
             * @return StructArray - an Array with the appropriate dimensions and fieldNames
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw DuplicateFieldNameInStructArrayException - if duplicate field names are specified
             */
            StructArray createStructArray(ArrayDimensions dims, std::vector<std::string> fieldNames) {
                NameList field_names(create_names(fieldNames.size()));
                for (const auto& field : fieldNames) {
                    add_name(field_names.getImpl(), field.c_str(), field.size());
                }
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_struct_array(pImpl.get(), &dims[0], dims.size(), field_names.getImpl(), &impl));
                return matlab::data::detail::Access::createObj<StructArray>(impl);
            }

            /**
             * Creates an Enumeration array. The class of the EnumArray is indicated by className.
             * The array is initialized with the list of enumeration names provided by the caller.
             *
             * @param dims array dimensions
             * @param className class name of the enumeration array
             * @param enums List of the enumeration names
             *
             * @return EnumArray - the EnumArray that was created
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw MustSpecifyClassNameException - if an empty class name is specified
             * @throw WrongNumberOfEnumsSuppliedException - if the wrong number of enumerations is provided
             */
            EnumArray createEnumArray(ArrayDimensions dims, std::string className, std::vector<std::string> enums) {
                NameList enum_names(create_names(enums.size()));
                for (const auto& e : enums) {
                    add_name(enum_names.getImpl(), e.c_str(), e.size());
                }
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_enum_array_with_enums(pImpl.get(), &dims[0], dims.size(), className.c_str(), className.size(), enum_names.getImpl(), &impl));
                return matlab::data::detail::Access::createObj<EnumArray>(impl);
            }

            /**
             * Creates an Enumeration array. The class of the EnumArray is indicated by className.
             *
             * @param dims array dimensions
             * @param className class name of the enumeration array
             *
             * @return EnumArray - the EnumArray that was created
             *
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             * @throw MustSpecifyClassNameException - if an empty class name is specified
             */
            EnumArray createEnumArray(ArrayDimensions dims, std::string className) {
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_enum_array(pImpl.get(), &dims[0], dims.size(), className.c_str(), className.size(), &impl));
                return matlab::data::detail::Access::createObj<EnumArray>(impl);
            }

            /**
             * Creates an empty Array. It contains no elements.
             *
             * @return Array - an empty Array
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */ 
            Array createEmptyArray() {
                return createArray<double>({0,0});
            }

            /**
             * Creates a buffer which can be passed into createArrayFromBuffer()
             * No data copies are made when creating an Array from a buffer. Data
             * must be in column major order.
             *
             * @param numberOfElements - the number of elements, not the actual buffer size
             *
             * @return buffer_ptr<T> - unique_ptr containing the buffer
             *
             * @throw InvalidArrayTypeException - if buffer type is not recognized as valid
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */ 
            template <typename T>
            buffer_ptr_t<T> createBuffer(size_t numberOfElements) {
                void* buffer = nullptr;
                buffer_deleter_t deleter = nullptr;                                                                      
                detail::throwIfError(create_buffer(pImpl.get(), &buffer, &deleter, static_cast<int>(GetArrayType<T>::type), numberOfElements));
                return buffer_ptr_t<T>(static_cast<T*>(buffer), deleter);
            }

            /**
             * Creates a TypedArray<T> using the given buffer
             * No data copies are made when creating an Array from a buffer. The TypedArray<T>
             * takes ownership of the buffer. Supplied data must be in column major order.
             *
             * @param dims - the dimensions for the Array
             * @param buffer - buffer containing the data. Buffer will not be copied. Array takes ownership
             *
             * @return TypedArray<T> - TypedArray<T> which wraps the buffer
             *
             * @throw InvalidArrayTypeException - if data type is not recognized as valid
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             */ 
            template<typename T>
            TypedArray<T> createArrayFromBuffer(matlab::data::ArrayDimensions dims,
                                                buffer_ptr_t<T> buffer) {
                buffer_deleter_t deleter = buffer.get_deleter();
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_array_from_buffer(pImpl.get(), 
                                                              static_cast<int>(GetArrayType<T>::type), 
                                                              &dims[0], 
                                                              dims.size(),
                                                              buffer.release(),
                                                              deleter,
                                                              &impl));

                return matlab::data::detail::Access::createObj<TypedArray<T>>(impl);
            }
            
            /**
             * Creates a SparseArray<T> with the given dimensions.
             * Only two dimensions are allowed for sparse arrays.
             *
             * @param dims - the dimensions for the Array
             * @param nnz - number of non-zero elements
             * @param data - buffer containing the non-zero elements. Buffer will not be copied. Array takes ownership
             * @param rows - buffer containing the row value for each element. Buffer will not be copied. Array takes ownership
             * @param cols - buffer containing the row value for each element. Buffer will not be copied. Array takes ownership
             *
             * @return SparseArray<T> - a sparse array with the appropriate type and dimensions, and data
             *
             * @throw InvalidArrayTypeException - if array type is not recognized as valid
             * @throw InvalidDimensionsInSparseArrayException if more than two dimensions are specified.
             * @throw matlab::OutOfMemoryException - if the array could not be allocated
             *
             */
            template<typename T>
            typename std::enable_if<std::is_same<double, T>::value ||
                                    std::is_same<bool, T>::value ||
                                    std::is_same<std::complex<double>, T>::value, SparseArray<T>>::type
            createSparseArray(matlab::data::ArrayDimensions dims, 
                              size_t nnz,
                              buffer_ptr_t<T> data, 
                              buffer_ptr_t<size_t> rows, 
                              buffer_ptr_t<size_t> cols) {

                // If more than two dimensions are specified, throw an exception
                if (dims.size() > 2) {
                    throw InvalidDimensionsInSparseArrayException(std::string("Sparse Array can only have 2 dimensions"));
                }
                
                buffer_deleter_t data_deleter = data.get_deleter();
                buffer_deleter_t rows_deleter = rows.get_deleter();
                buffer_deleter_t cols_deleter = cols.get_deleter();
                
                matlab::data::impl::ArrayImpl* impl = nullptr;
                detail::throwIfError(create_sparse_array_from_buffer(pImpl.get(), 
                                                                     static_cast<int>(GetSparseArrayType<T>::type), 
                                                                     &dims[0], 
                                                                     dims.size(),
                                                                     nnz,
                                                                     data.release(),
                                                                     data_deleter,
                                                                     rows.release(),
                                                                     rows_deleter,
                                                                     cols.release(),
                                                                     cols_deleter,
                                                                     &impl));
                return matlab::data::detail::Access::createObj<SparseArray<T>>(impl);
            }
       
          protected:

            std::shared_ptr<matlab::data::impl::ArrayFactoryImpl> pImpl;

          private:
            ArrayFactory& operator=(ArrayFactory const& rhs) = delete;
            ArrayFactory(const ArrayFactory &rhs) = delete;

            class NameList {
              public:
                NameList(detail::NameListImpl* impl) :
                    pImpl(std::shared_ptr<detail::NameListImpl>
                          (impl, [](detail::NameListImpl* ptr) {
                              names_destroy_impl(ptr);
                          })) {}

                detail::NameListImpl* getImpl() const {
                    return pImpl.get();
                }

              private:
                std::shared_ptr<detail::NameListImpl> pImpl;
            };
            template<typename T> 
            typename std::enable_if<!std::is_base_of<Array, T>::value, Array>::type createValue(T value) {
                return createScalar(std::move(value));
            }
            template<typename T> 
            typename std::enable_if<std::is_base_of<Array, T>::value, Array>::type createValue(T value) {
                return value;
            }
            template <typename T>
            void setCellValues(TypedIterator<Array> it, TypedIterator<Array> end, T value) {
                *it = createValue(std::move(value));
            }
            template <typename T, typename... Targs>
            void setCellValues(TypedIterator<Array> it, TypedIterator<Array> end, T value, Targs... args) {
                *it = createValue(std::move(value));
                if (++it != end) {
                    setCellValues(std::move(it), std::move(end), args...);
                }
            }
        };
    }
}

#endif
            
            
