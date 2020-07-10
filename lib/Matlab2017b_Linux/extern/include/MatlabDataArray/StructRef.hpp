/* Copyright 2015-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_STRUCT_REF_HPP
#define MATLAB_EXTDATA_STRUCT_REF_HPP

#include "TypedArrayRef.hpp"
#include "Struct.hpp"
#include "MDArray.hpp"
#include "Exception.hpp"

#include "detail/struct_interface.hpp"
#include "detail/array_interface.hpp"
#include "detail/HelperFunctions.hpp"
#include "detail/ReferenceHelpers.hpp"

namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;
            class Access;
        }

        
        template<>
        class Reference<Struct> : public ReferenceHolder {
          public:
            
            typedef TypedIterator<Array> iterator;
            typedef TypedIterator<Array> const_iterator;
            /**
             * operator[] index into the Struct with a field name
             *
             * @param idx - the field name
             * @return Reference<Array> - reference to the field value
             * @throw InvalidFieldNameException if the field does not exist in the struct
             *
             */
            Reference<Array> operator[](std::string idx) {
                return detail::Access::createObj<Reference<Array>>(detail::getRef(pImpl.get(),
                                                                                  idx.c_str(),
                                                                                  idx.size(),
                                                                                  true));
            }

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
             * Return an iterator to the beginning of the list of fields
             *
             * @return iterator
             * @throw none
             */
            iterator begin() MW_NOEXCEPT {
                return detail::Access::createObj<iterator>(array_reference_begin(pImpl.get(), true));
            }

            /**
             * Return an iterator to the end of the list of fields
             *
             * @return iterator
             * @throw none
             */
            iterator end() MW_NOEXCEPT {
                return detail::Access::createObj<iterator>(array_reference_end(pImpl.get(), true));
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

            /**
             * cast to a Struct
             *
             * @return Struct
             * @throw none
             */
            operator Struct() const MW_NOEXCEPT {
                return detail::Access::createObj<Struct>(pImpl);
            }
            
          private:

            friend class detail::Access;

            Reference(detail::ReferenceImpl* impl) MW_NOEXCEPT :
                ReferenceHolder(impl) {}

            Reference() = delete;

            Reference(std::shared_ptr<detail::ReferenceImpl> impl) MW_NOEXCEPT :
                ReferenceHolder(std::move(impl)) {}

        };
    }
}

#endif
