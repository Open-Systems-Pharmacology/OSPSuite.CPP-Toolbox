/* Copyright 2014-2017 The MathWorks, Inc. */

#ifndef MATLAB_DATA_STRUCT_HPP
#define MATLAB_DATA_STRUCT_HPP

#include "MDArray.hpp"
#include "TypedIterator.hpp"

#include "detail/struct_interface.hpp"
#include "detail/reference_interface.hpp"
#include "detail/ExceptionHelpers.hpp"
#include "detail/HelperFunctions.hpp"

namespace matlab {
    namespace data {        
        namespace detail {
            class ReferenceImpl;
            class Access;
        }
        namespace impl {
            class ArrayImpl;
        }
        /**
         * Struct Element - ElementType for StructArrays. Can be used to inspect the fields of the Struct, but not
         * modify the Struct.
         */
        class Struct {
          public:
            typedef TypedIterator<Array const> const_iterator;

            /**
             * return a shared copy of the Array found at the field specified
             *
             * @param idx - field name
             * @return Array - shared copy of the Array found at the field specified
             * @throw InvalidFieldNameException if the field does not exist in the struct
             */
            Array operator[](std::string idx) const {
                matlab::data::impl::ArrayImpl* arrayImpl = nullptr;
                detail::throwIfError(struct_get_field_copy(pImpl.get(), idx.c_str(), idx.size(), &arrayImpl));
                return detail::Access::createObj<Array>(arrayImpl);
            }
            
            /**
             * Return a const_iterator to the beginning of the list of fields
             *
             * @return const_iterator
             * @throw none
             */
            const_iterator begin() const MW_NOEXCEPT {
                return cbegin();
            }

            /**
             * Return a const_iterator to the beginning of the list of fields
             *
             * @return const_iterator
             * @throw none
             */
            const_iterator cbegin() const MW_NOEXCEPT {
                return detail::Access::createObj<const_iterator>(array_reference_begin(pImpl.get(), false));
            }

            /**
             * Return a const_iterator to the end of the list of fields
             *
             * @return const_iterator
             * @throw none
             */
            const_iterator end() const MW_NOEXCEPT {
                return cend();
            }

            /**
             * Return a const_iterator to the end of the list of fields
             *
             * @return const_iterator
             * @throw none
             */
            const_iterator cend() const MW_NOEXCEPT {
                return detail::Access::createObj<const_iterator>(array_reference_end(pImpl.get(), false));
            }

          private:
            Struct() = delete;

            Struct(std::shared_ptr<detail::ReferenceImpl> impl) MW_NOEXCEPT :
                pImpl(impl) {}

            std::shared_ptr<detail::ReferenceImpl> pImpl;
            friend class detail::Access;
        };
    }
}

#endif
