/* Copyright 2014-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_ARRAY_REF_EXT_HPP
#define MATLAB_EXTDATA_ARRAY_REF_EXT_HPP

#include "MDArray.hpp"
#include "Enumeration.hpp"
#include "String.hpp"
#include "ReferenceHolder.hpp"

#include "detail/reference_interface.hpp"
#include "detail/HelperFunctions.hpp"
#include "detail/ReferenceHelpers.hpp"
#include "detail/ExceptionHelpers.hpp"

#include <memory>
#include <type_traits>

namespace matlab {
    namespace data {
        namespace detail {
            class Access;
        }
        /**
         * The ArrayReferenceExt class is a base class for all reference types that refer to Arrays.
         * It provides basic Array information.
         */
        class ArrayReferenceExt : public ReferenceHolder {
          public:
            /**
             * Returns the ArrayType of the array
             *
             * @return ArrayType the ArrayType
             *
             * @throw InvalidArrayIndexException - if one of the indices is out of range for this Array
             * @throw NotEnoughIndicesProvidedException - if not enough indices were provided
             * @throw InvalidArrayTypeException - if data type returned is not recognized
             */
            ArrayType getType() const {
                int retVal;
                detail::throwIfError(array_reference_get_type(pImpl.get(), 
                                                              &retVal));
                return ArrayType(retVal);
            }
            
            /**
             * Get the array's dimensions.
             *
             * @return ArrayDimensions
             *
             * @throw InvalidArrayIndexException - if one of the indices is out of range for this Array
             * @throw NotEnoughIndicesProvidedException - if not enough indices were provided
             */
            ArrayDimensions getDimensions() const {
                size_t numDims = 0;
                size_t* dims = nullptr;
                detail::throwIfError(array_reference_get_dimensions(pImpl.get(), &numDims, &dims));
                return ArrayDimensions(dims, dims+numDims);
            }

            /**
             * Get the number of elements in the Array that this reference refers to
             *
             * @return the number of elements in the Array
             *
             * @throw InvalidArrayIndexException - if one of the indices is out of range for this Array
             * @throw NotEnoughIndicesProvidedException - if not enough indices were provided
             */
            size_t getNumberOfElements() const {
                size_t retVal;
                detail::throwIfError(array_reference_get_num_elements(pImpl.get(), &retVal));
                return retVal;
            }

            /**
             * Determine if this array is empty
             *
             * @return bool - true if empty, otherwise false
             *
             * @throw InvalidArrayIndexException - if one of the indices is out of range for this Array
             * @throw NotEnoughIndicesProvidedException - if not enough indices were provided
             */
            bool isEmpty() const {
                bool retVal;
                detail::throwIfError(array_reference_is_empty(pImpl.get(), &retVal));
                return retVal;
            }

          protected:
            friend class detail::Access;
            
            ArrayReferenceExt(detail::ReferenceImpl* impl) MW_NOEXCEPT :
                ReferenceHolder(impl) {}

            ArrayReferenceExt(std::shared_ptr<detail::ReferenceImpl> impl) MW_NOEXCEPT :
                ReferenceHolder(impl) {}

            ArrayReferenceExt(ArrayReferenceExt&& rhs) MW_NOEXCEPT :
                ReferenceHolder(std::move(rhs)) {}

            ArrayReferenceExt& operator= (ArrayReferenceExt&& rhs) MW_NOEXCEPT {
                ReferenceHolder::operator= (std::move(rhs));
                return *this;
            }

            ArrayReferenceExt(const ArrayReferenceExt& rhs) MW_NOEXCEPT :
                ReferenceHolder(rhs) {}

            ArrayReferenceExt& operator= (const ArrayReferenceExt& rhs) MW_NOEXCEPT {
                ReferenceHolder::operator= (rhs);
                return *this;
            }

          private:
            ArrayReferenceExt() = delete;

        };

        /**
         * The IndexableArrayRef class is a base class for Reference<Array> to enable the reference
         * to be further indexed. 
         */
        template <bool is_const_ref>
        class IndexableArrayRef : public ArrayReferenceExt {
          public:

            /**
             * Enables [] indexing on a reference to an Array.
             *
             * The return value ArrayElementRef<is_const_ref> allows the element of the array to be
             * modified or retrieved: For example:
             *     auto cell_arr = factory.createArray<Array>({2,2});
             *     cell_arr[0][0] = factory.createArray<double>({3,3});
             *     Reference<Array> ref = cell_arr[0][0];
             *     ref[0][0] = 5.5;
             *
             * @param idx - the first array index
             *
             * @return ArrayElementRef<is_const_ref> - contains the index specified
             *
             * @throw InvalidFieldNameException if the Reference<Array> was created with an invalid fieldname for a struct
             */
            ArrayElementRef<is_const_ref> operator[](size_t idx) {
                detail::ReferenceImpl* newRef = nullptr;
                detail::throwIfError(array_reference_get_reference(pImpl.get(), !is_const_ref, idx, &newRef));
                return detail::Access::createObj<ArrayElementRef<is_const_ref>>(newRef);
            }

            /**
             * Enables [] indexing on a reference to an Array.
             *
             * The return value ArrayElementRef<true> allows the element of the array to be
             * modified or retrieved: For example:
             *     auto cell_arr = factory.createArray<Array>({2,2});
             *     cell_arr[0][0] = factory.createScalar(5.5);
             *     Reference<Array const> ref = cell_arr[0][0];
             *     double val = ref[0];
             *
             * @param idx - the first array index
             *
             * @return ArrayElementRef<true> - contains the index specified
             *
             * @throw InvalidFieldNameException if the Reference<Array> was created with an invalid fieldname for a struct
             */
            ArrayElementRef<true> operator[](size_t idx) const {
                detail::ReferenceImpl* newRef = nullptr;
                detail::throwIfError(array_reference_get_reference(pImpl.get(), false, idx, &newRef));
                return detail::Access::createObj<ArrayElementRef<true>>(newRef);
            }

          protected:
            IndexableArrayRef(detail::ReferenceImpl* impl) MW_NOEXCEPT :
                ArrayReferenceExt(impl) {}

            IndexableArrayRef(std::shared_ptr<detail::ReferenceImpl> impl) MW_NOEXCEPT :
                ArrayReferenceExt(impl) {}

            IndexableArrayRef(IndexableArrayRef<is_const_ref>&& rhs) MW_NOEXCEPT :
                ArrayReferenceExt(std::move(rhs)) {}

            IndexableArrayRef<is_const_ref>& operator= (IndexableArrayRef<is_const_ref>&& rhs) MW_NOEXCEPT {
                ArrayReferenceExt::operator= (std::move(rhs));
                return *this;
            }

            IndexableArrayRef(const IndexableArrayRef<is_const_ref>& rhs) MW_NOEXCEPT :
                ArrayReferenceExt(rhs) {}

            IndexableArrayRef<is_const_ref>& operator= (const IndexableArrayRef<is_const_ref>& rhs) MW_NOEXCEPT {
                ArrayReferenceExt::operator= (rhs);
                return *this;
            }
        };
    }
}

#endif

