/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_REFERENCE_HELPERS_HPP
#define MATLAB_EXTDATA_REFERENCE_HELPERS_HPP

#include "reference_interface.hpp"
#include "string_interface.hpp"
#include "struct_interface.hpp"
#include "exception_interface.hpp"
#include "ExceptionHelpers.hpp"

#include <stdexcept>

namespace matlab {
    namespace data {
        template <typename T> class Reference;
        
        namespace detail {

            template <typename T, typename U>
            inline typename std::enable_if<std::is_base_of<Array, U>::value>::type validateReference(ReferenceImpl* impl, bool unshare) {
                detail::throwIfError(typed_array_reference_validate_type(impl, static_cast<int>(U::type), false));
            }
                
            template <typename T, typename U>
            inline typename std::enable_if<!std::is_base_of<Array, U>::value>::type validateReference(ReferenceImpl* impl, bool unshare) {
                static_assert(std::is_same<T, U>::value, "Reference type must match array type");
                detail::throwIfError(reference_validate_index(impl));
            }

            template <typename T>
            inline typename std::enable_if<std::is_base_of<Array, T>::value>::type validateUntypedReference(ReferenceImpl* impl) {
                detail::throwIfError(typed_array_reference_validate_type(impl, static_cast<int>(T::type), false));
            }
                
            template <typename T>
            inline typename std::enable_if<!std::is_base_of<Array, T>::value>::type validateUntypedReference(ReferenceImpl* impl) {
                detail::throwIfError(reference_validate_type(impl, static_cast<int>(GetArrayType<T>::type)));
            }

            template <typename T>
            inline typename std::enable_if<!std::is_base_of<Array, T>::value>::type validateUntypedArrayReference(ReferenceImpl* impl) {
                detail::throwIfError(typed_array_reference_validate_type(impl, static_cast<int>(GetArrayType<T>::type), false));
            }
            template <typename T>
            inline typename std::enable_if<std::is_base_of<Array, T>::value>::type validateUntypedArrayReference(ReferenceImpl* impl) {
                detail::throwIfError(typed_array_reference_validate_type(impl, static_cast<int>(GetArrayType<T>::type), false));
            }

            template <typename T>
            inline T castToType(void* value, ArrayType type) {
                switch (type) {
                  case ArrayType::LOGICAL: return static_cast<T>(*static_cast<bool*>(value));
                  case ArrayType::DOUBLE: return static_cast<T>(*static_cast<double*>(value));
                  case ArrayType::SINGLE: return static_cast<T>(*static_cast<float*>(value));
                  case ArrayType::INT64: return static_cast<T>(*static_cast<int64_t*>(value));
                  case ArrayType::INT32: return static_cast<T>(*static_cast<int32_t*>(value));
                  case ArrayType::INT16: return static_cast<T>(*static_cast<int16_t*>(value));
                  case ArrayType::INT8: return static_cast<T>(*static_cast<int8_t*>(value));
                  case ArrayType::UINT64: return static_cast<T>(*static_cast<uint64_t*>(value));
                  case ArrayType::UINT32: return static_cast<T>(*static_cast<uint32_t*>(value));
                  case ArrayType::UINT16: return static_cast<T>(*static_cast<uint16_t*>(value));
                  case ArrayType::UINT8: return static_cast<T>(*static_cast<uint8_t*>(value));
                  default:
                    throw TypeMismatchException(std::string("Can't convert this element"));
                }
                return T();
            }

            template <>
            inline bool castToType<bool>(void* value, ArrayType type) {
                switch (type) {
                  case ArrayType::LOGICAL: return *static_cast<bool*>(value);
                  case ArrayType::DOUBLE: return *static_cast<double*>(value) != 0;
                  case ArrayType::SINGLE: return *static_cast<float*>(value) != 0;
                  case ArrayType::INT64: return *static_cast<int64_t*>(value) != 0;
                  case ArrayType::INT32: return *static_cast<int32_t*>(value) != 0;
                  case ArrayType::INT16: return *static_cast<int16_t*>(value) != 0;
                  case ArrayType::INT8: return *static_cast<int8_t*>(value) != 0;
                  case ArrayType::UINT64: return *static_cast<uint64_t*>(value) != 0;
                  case ArrayType::UINT32: return *static_cast<uint32_t*>(value) != 0;
                  case ArrayType::UINT16: return *static_cast<uint16_t*>(value) != 0;
                  case ArrayType::UINT8: return *static_cast<uint8_t*>(value) != 0;
                  default:
                    throw TypeMismatchException(std::string("Can't convert this element"));
                }
                return false;
            }

            template <>
            inline CHAR16_T castToType<CHAR16_T>(void* value, ArrayType type) {
                switch (type) {
                  case ArrayType::CHAR: return *static_cast<CHAR16_T*>(value);
                  default:
                    throw TypeMismatchException(std::string("Can't convert this element"));
                }
                return false;
            }

            template <typename U>
            inline typename std::enable_if<std::is_arithmetic<U>::value, U&>::type getElement(std::shared_ptr<ReferenceImpl> impl) {
                void* value = nullptr;
                throwIfError(typed_reference_get_pod_value(impl.get(),
                                                           &value));
                return *(static_cast<U*>(value));
            }

            template <typename U>
            inline typename std::enable_if<is_complex<U>::value, U>::type getElement(std::shared_ptr<ReferenceImpl> impl) {                
                void* real = nullptr;
                void* imag = nullptr;
                throwIfError(typed_reference_get_complex_value(impl.get(),
                                                               &real,
                                                               &imag));

                typename U::value_type r = *(static_cast<typename U::value_type*>(real));
                typename U::value_type i = *(static_cast<typename U::value_type*>(imag));

                return U(r,i);
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value, MATLABString>::type getElement(std::shared_ptr<ReferenceImpl> impl) {
                char16_t* str = nullptr;
                size_t strlen = 0;
                throwIfError(string_get_value(impl.get(), &str, &strlen));
                if (str != nullptr) {
                    return MATLABString(String(str, strlen));
                }

                return MATLABString();
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<Array, U>::value, Array>::type getElement(std::shared_ptr<ReferenceImpl> impl) {
                impl::ArrayImpl* arr_impl = nullptr;
                array_reference_shared_copy(impl.get(), &arr_impl);
                return detail::Access::createObj<U>(arr_impl);
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<Struct, U>::value, Struct>::type getElement(std::shared_ptr<ReferenceImpl> impl) {
                return detail::Access::createObj<U>(impl);
            }

            template <typename U>
            inline typename std::enable_if<std::is_arithmetic<U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                throwIfError(typed_reference_set_pod_value(impl,
                                                           type,
                                                           &rhs));
            }

            template <typename U>
            inline typename std::enable_if<is_complex<U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                typename U::value_type real = rhs.real();
                typename U::value_type imag = rhs.imag();

                throwIfError(typed_reference_set_complex_value(impl,
                                                               type,
                                                               &real,
                                                               &imag));
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<String, U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                throwIfError(reference_set_char16_string(impl, rhs.c_str(), rhs.size()));
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<std::string, U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                throwIfError(reference_set_string(impl, rhs.c_str(), rhs.size()));
            }

            template <typename U>
            inline typename std::enable_if<std::is_base_of<Array, U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                throwIfError(reference_set_reference_value(impl,
                                                           detail::Access::getImpl<matlab::data::impl::ArrayImpl>(rhs)));
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value>::type setElement(ReferenceImpl* impl, U rhs, int type) {
                if (rhs) {
                    String elem = *rhs;
                    throwIfError(reference_set_char16_string(impl, elem.c_str(), elem.size()));
                }
                else {
                    throwIfError(reference_set_missing_char16_string(impl));
                }
            }

            template <typename U>
            inline typename std::enable_if<!std::is_same<Enumeration, U>::value &&
                                           !std::is_same<MATLABString, U>::value, std::string>::type getString(ReferenceImpl* impl) {
                throwIfError(static_cast<int>(ExceptionType::InvalidDataType));
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<Enumeration, U>::value, std::string>::type getString(ReferenceImpl* impl) {
                char* str = nullptr;
                size_t strlen = 0;
                throwIfError(enum_get_value(impl,
                                            &str,
                                            &strlen));
                return std::string(str, strlen);
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value, std::string>::type getString(ReferenceImpl* impl) {
                char16_t* str = nullptr;
                size_t strlen = 0;
                bool missing = false;
                throwIfError(string_get_value(impl,
                    &str,
                    &strlen));

                if (str == nullptr) {
                    throw std::runtime_error("Missing string");
                }
                String temp(str, strlen);
                if (!detail::isAscii7(temp)) {
                    throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                }
                return std::string(temp.begin(), temp.end());
            }

            template <typename U>
            inline typename std::enable_if<!std::is_same<Enumeration, U>::value &&
                                           !std::is_same<MATLABString, U>::value>::type setString(ReferenceImpl* impl, std::string rhs) {
                throwIfError(static_cast<int>(ExceptionType::InvalidDataType));
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<Enumeration, U>::value>::type setString(ReferenceImpl* impl, std::string rhs) {
                throwIfError(enum_set_value(impl, 
                                            const_cast<char*>(rhs.c_str()),
                                            rhs.size()));
            }
            
            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value>::type setString(ReferenceImpl* impl, std::string rhs) {
                if (!detail::isAscii7(rhs)) {
                    throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                }
                String tmp(rhs.begin(), rhs.end());
                setElement<String>(impl, std::move(tmp), 0);
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value>::type setString(ReferenceImpl* impl, String rhs) {
                setElement<String>(impl, std::move(rhs), 0);
            }

            template <typename U>
            inline typename std::enable_if<std::is_same<MATLABString, U>::value>::type setString(ReferenceImpl* impl, MATLABString rhs) {
                setElement<MATLABString>(impl, std::move(rhs), 0);
            }

            inline ReferenceImpl* getRef(ReferenceImpl* impl, const char* field, size_t len, bool unshare) {
                ReferenceImpl* retVal = nullptr;
                size_t numIdx = 0;
                throwIfError(struct_reference_get_index(impl, field, len, &numIdx));
                reference_get_reference_value(impl, unshare, &retVal);
                reference_add_index(retVal, numIdx);
                return retVal;
            }

            template <typename T, bool is_const_ref>
            typename std::enable_if<std::is_arithmetic<T>::value, T>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                void* value = nullptr;
                detail::throwIfError(reference_get_pod_value(impl.get(),
                                                             &value));
                int type;
                detail::throwIfError(reference_get_container_type(impl.get(), &type));

                return detail::castToType<T>(value, ArrayType(type));
            }

            template <typename T, bool is_const_ref>
            typename std::enable_if<std::is_base_of<Array, T>::value, T>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                detail::throwIfError(typed_array_reference_validate_type(impl.get(), static_cast<int>(T::type), false));
                impl::ArrayImpl* arr_impl = nullptr;
                array_reference_shared_copy(impl.get(), &arr_impl);
                return detail::Access::createObj<T>(arr_impl);
            }


            template <typename T, bool is_const_ref>
            typename std::enable_if<is_complex<T>::value, T>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                void* real = nullptr;
                void* imag = nullptr;
                detail::throwIfError(reference_get_complex_value(impl.get(),
                                                                 &real,
                                                                 &imag));

                typename T::value_type r = *(static_cast<typename T::value_type*>(real));
                typename T::value_type i = *(static_cast<typename T::value_type*>(imag));

                return T(r,i);
            }

            template <typename T, bool is_const_ref>
            typename std::enable_if<std::is_same<T, String>::value, String>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                int type;
                detail::throwIfError(reference_get_container_type(impl.get(), &type));
                if (ArrayType(type) == ArrayType::MATLAB_STRING) {
                    char16_t* str = nullptr;
                    size_t strlen = 0;
                    detail::throwIfError(string_get_value(impl.get(),
                        &str,
                        &strlen));
                    return String(str, strlen);
                }
                throw TypeMismatchException(std::string("Can't convert this element to a matlab::data::String"));
            }

            template <typename T, bool is_const_ref>
            typename std::enable_if<std::is_base_of<ReferenceHolder, T>::value, T>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                static_assert(!is_const_ref, "Can only get a reference to a non const ref");
                validateUntypedReference<typename T::ref_type>(impl.get());
                return detail::Access::createObj<T>(impl);
            }

            template <typename T, bool is_const_ref>
            typename std::enable_if<std::is_same<T, std::string>::value, std::string>::type castTo(std::shared_ptr<detail::ReferenceImpl> impl) {
                int arrType;
                detail::throwIfError(reference_get_container_type(impl.get(), &arrType));

                if (ArrayType(arrType) == ArrayType::ENUM) {
                    char* str = nullptr;
                    size_t strlen = 0;
                    detail::throwIfError(enum_get_value(impl.get(),
                                                        &str,
                                                        &strlen));
                    return std::string(str, strlen);
                }
                else if (ArrayType(arrType) == ArrayType::MATLAB_STRING) {
                    char16_t* str = nullptr;
                    size_t strlen = 0;
                    detail::throwIfError(string_get_value(impl.get(),
                        &str,
                        &strlen));
                    String temp(str, strlen);
                    if (!detail::isAscii7(temp)) {
                        throw NonAsciiCharInInputDataException(std::string("Input data can only contain ASCII characters"));
                    }
                    return std::string(temp.begin(), temp.end());
                }
                throw TypeMismatchException(std::string("Can't convert this element to a std::string"));
            }
        }
    }
}

#endif
