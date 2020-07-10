/* Copyright 2014 - 2016 The MathWorks, Inc. */

#ifndef MATLAB_DATA_FIELDIDENTIFIER_HPP
#define MATLAB_DATA_FIELDIDENTIFIER_HPP


#include "detail/field_id_interface.hpp"

#include <string>
#include <memory>

namespace matlab {
    namespace data {
        namespace detail {
            class RefCounted;
            class Access;
        }

        template <typename T> class ForwardIterator;

        /**
         * The MATLABFieldIdentifier is a class used to indentify field names
         */
        class MATLABFieldIdentifier {
          public:

            /**
             * Construct a MATLABFieldIdentifier
             *
             * @return newly constructed MATLABFieldIdentifer
             * @throw none
             */
            MATLABFieldIdentifier() MW_NOEXCEPT {}

            /**
             * Construct a MATLABFieldIdentifier from a std::string
             * @param str - the string that contains the field name
             * @return newly constructed MATLABFieldIdentifer
             * @throw none
             */
            MATLABFieldIdentifier(std::string str) MW_NOEXCEPT :
                pImpl(std::shared_ptr<matlab::data::detail::RefCounted>
                      (field_id_create_impl(str.c_str(), str.size()), [](matlab::data::detail::RefCounted* ptr) {
                          field_id_destroy_impl(ptr);
                      })) {}
                      
            /**
             * Move assignment operator
             * @param rhs - the MATLABFieldIdentifier to be moved
             * @return updated MATLABFieldIdentifer
             * @throw none
             */
            MATLABFieldIdentifier& operator=(MATLABFieldIdentifier&& rhs) MW_NOEXCEPT {
                pImpl = std::move(rhs.pImpl);
                return *this;
            }

            /**
             * Move constructor
             * @param rhs - the MATLABFieldIdentifier to be moved
             * @return newly constructed MATLABFieldIdentifer
             * @throw none
             */
            MATLABFieldIdentifier(MATLABFieldIdentifier&& rhs) MW_NOEXCEPT :
                pImpl(std::move(rhs.pImpl)) {}

            /**
             * Assignment operator
             * @param rhs - the MATLABFieldIdentifier to be copied
             * @return updated MATLABFieldIdentifier which is a shared copy of the original MATLABFieldIdentifier
             * @throw none
             */
            MATLABFieldIdentifier& operator=(MATLABFieldIdentifier const& rhs) MW_NOEXCEPT {
                pImpl = rhs.pImpl;
                return *this;
            }

            /**
             * copy constructor
             * @param rhs - the MATLABFieldIdentifier to be copied
             * @return newly constructed MATLABFieldIdentifer
             * @throw none
             */
            MATLABFieldIdentifier(const MATLABFieldIdentifier &rhs) MW_NOEXCEPT :
                pImpl(rhs.pImpl) {}

            /**
             * Check to see if the 2 MATLABFieldIdentifiers are identical
             * @param rhs - MATLABFieldIdentifier to be compared
             * @return true if the MATLABFieldIdentifiers are identical
             * @throw none
             */
            bool operator==( const MATLABFieldIdentifier& rhs) const MW_NOEXCEPT {
                return field_id_is_equal(pImpl.get(), rhs.pImpl.get());
            }

            /**
             * Get a std::string representation of the MATLABFieldIdentifier
             * @return std::string - std::string representation of the MATLABFieldIdentifier
             * @throw none
             */
            operator std::string() const MW_NOEXCEPT {
                const char* str = nullptr;
                size_t len = 0;
                field_id_get_string(pImpl.get(), &str, &len);
                return std::string(str, len);
            }

          private:
            std::shared_ptr<matlab::data::detail::RefCounted> pImpl;

            friend class detail::Access;

            MATLABFieldIdentifier(matlab::data::detail::RefCounted* impl) MW_NOEXCEPT :
                pImpl(std::shared_ptr<matlab::data::detail::RefCounted>
                      (impl, [](matlab::data::detail::RefCounted* ptr) {
                          field_id_destroy_impl(ptr);
                      })) {}
        };
    }
}
#endif
