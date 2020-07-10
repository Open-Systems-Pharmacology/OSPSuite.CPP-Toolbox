/* Copyright 2014-2017 The MathWorks, Inc. */

#ifndef MDA_ARRAY_HPP_
#define MDA_ARRAY_HPP_


#include "matlab_data_array_defs.hpp"
#include "matlab_extdata_defs.hpp"
#include "ArrayDimensions.hpp"
#include "ArrayType.hpp"
#include "ArrayElementRef.hpp"

#include "detail/array_interface.hpp"
#include "detail/HelperFunctions.hpp"

#include <functional>
#include <stdint.h>
#include <exception>


namespace matlab {
    namespace data {

        template<template<class> class IteratorType, class ElementType> class Range;
        template<typename T> class TypedIterator;
        
        namespace impl {
            class ArrayImpl;
        }
        namespace detail {
            class Access;
        }
        
        /**
         * Array provides an API to access general information about an
         * array.  Specific subclasses will provide access to type-specific
         * data for different kinds of arrays.  
         */
        class Array {
          public:

            static const ArrayType type = ArrayType::UNKNOWN;
            
            /**
             * Default constructor - creates an empty Array
             * @return - newly constructed Array
             * @throw - none
             */
            Array() MW_NOEXCEPT :
                pImpl(std::shared_ptr<matlab::data::impl::ArrayImpl>
                      (array_create_empty(), [](matlab::data::impl::ArrayImpl* ptr) {
                          array_destroy_impl(ptr);
                      })) {}

            /**
             * Destructor
             *
             * @throw none
             */
            virtual ~Array() MW_NOEXCEPT
            {}
            
            /**
             * Move constructor
             *
             * @param - rhs Array value to be moved
             * @return - newly constructed Array
             * @throw none
             */
            Array(Array&& rhs) MW_NOEXCEPT :
                pImpl(std::move(rhs.pImpl)) {}
            
            /**
             * Move assignment operator
             *
             * @param - rhs Array value to be moved
             * @return - the updated Array
             * @throw none
             */
            Array& operator=(Array&& rhs) MW_NOEXCEPT {
                pImpl = std::move(rhs.pImpl);
                return *this;
            }
            
            /**
             * Assignment operator. The updated Array becomes a shared copy of the input Array
             *
             * @param - rhs Array value to be copied
             * @return - the updated Array
             * @throw none
             */
            Array& operator=(const Array& rhs) MW_NOEXCEPT {
                pImpl = rhs.pImpl;
                return *this;
            }

            /**
             * Copy consructor. The newly constructed Array becomes a shared copy of the input Array
             *
             * @param - rhs Array value to be copied
             * @return - new constructed Array
             * @throw none
             */
            Array(const Array &rhs) MW_NOEXCEPT :
                pImpl(rhs.pImpl) {}

            /**
             * Returns the ArrayType of the array
             *
             * @return - the ArrayType
             * @throw InvalidArrayType - if array type returned is not recognized as valid
             */
            ArrayType getType() const {
                int type;
                detail::throwIfError(array_get_type(pImpl.get(), &type));
                return static_cast<ArrayType>(type);
            }

            /**
             * Get the array's dimensions.
             *
             * @return ArrayDimensions vector of each dimension.
             * @throw none
             */
            ArrayDimensions getDimensions() const MW_NOEXCEPT {
                size_t numDims = 0;
                size_t* dims = nullptr;
                array_get_dimensions(pImpl.get(), &numDims, &dims);
                return ArrayDimensions(dims, dims+numDims);
            }

            /**
             * Get the number of elements in this array
             *
             * @return the number of elements in the array
             * @throw none
             */
            size_t getNumberOfElements() const MW_NOEXCEPT {
                return array_get_number_of_elements(pImpl.get());
            }

            /**
             * Determine if this is an empty array
             *
             * @return true if the array is empty
             * @throw none
             */
            bool isEmpty() const MW_NOEXCEPT {
                return array_is_empty(pImpl.get());
            }

            /**
             * Enables [] indexing on a array.
             *
             * The return value ArrayElementRef<false> allows the element of the array to be
             * modified or retrieved: For example:
             *     arr[1][1] = 5.5;
             *     double val = arr[0][3];
             *
             * @param idx - the first array index
             * @return ArrayElementRef<false> - contains the index specified
             * @throw InvalidArrayIndexException - if the index is invalid
             */
            ArrayElementRef<false> operator[](size_t idx) {
                impl::ArrayImpl* newImpl = nullptr;
                if (array_unshare(pImpl.get(), pImpl.unique(), &newImpl)) {
                    pImpl.reset(newImpl, [](impl::ArrayImpl* ptr) {
                        array_destroy_impl(ptr);
                    });
                }
                matlab::data::detail::ReferenceImpl* impl = nullptr;
                detail::throwIfError(array_create_reference(pImpl.get(), idx, &impl));
                return detail::Access::createObj<ArrayElementRef<false>>(impl);
            }

            /**
             * Enables [] indexing on a const array.
             *
             * The return value ArrayElementRef<true> allows the element of the array to be
             * retrieved, but not modified: For example:
             *     double val = arr[0][3];
             *
             * @param idx - the first array index
             * @return ArrayElementRef<true> - contains the index specified
             * @throw InvalidArrayIndexException - if the index is invalid
             */
            ArrayElementRef<true> operator[](size_t idx) const {
                matlab::data::detail::ReferenceImpl* impl = nullptr;
                detail::throwIfError(array_create_reference(pImpl.get(), idx, &impl));
                return detail::Access::createObj<ArrayElementRef<true>>(impl);
            }

          protected:
            Array(matlab::data::impl::ArrayImpl* impl) MW_NOEXCEPT :
                pImpl(std::shared_ptr<matlab::data::impl::ArrayImpl>
                      (impl, [](matlab::data::impl::ArrayImpl* ptr) {
                          array_destroy_impl(ptr);
                      })) {}

            std::shared_ptr<matlab::data::impl::ArrayImpl> pImpl;
            
            friend class detail::Access;

            template <typename T>
            friend Range<TypedIterator, T> getWritableElements(Array& a);
        };
    }
}

#endif
