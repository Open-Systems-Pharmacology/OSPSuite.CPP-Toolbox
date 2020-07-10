/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_ARRAY_REFERENCE_HOLDER_HPP
#define MATLAB_EXTDATA_ARRAY_REFERENCE_HOLDER_HPP

#include "detail/reference_interface.hpp"

#include <memory>


namespace matlab {
    namespace data {
        namespace detail {
            class ReferenceImpl;
        }

        class ReferenceHolder {
          protected:
            std::shared_ptr<detail::ReferenceImpl> pImpl;

            ReferenceHolder(detail::ReferenceImpl* impl) MW_NOEXCEPT :
                pImpl(std::shared_ptr<detail::ReferenceImpl>
                      (impl, [](detail::ReferenceImpl* ptr) {
                          reference_destroy_impl(ptr);
                      })) {}

            ReferenceHolder(std::shared_ptr<detail::ReferenceImpl> impl) MW_NOEXCEPT :
                pImpl(impl) {}

            ReferenceHolder(ReferenceHolder&& rhs) MW_NOEXCEPT :
                pImpl(std::move(rhs.pImpl)) {}

            ReferenceHolder& operator= (ReferenceHolder&& rhs) MW_NOEXCEPT {                            
                pImpl = std::move(rhs.pImpl);
                return *this;
            }
            
            ReferenceHolder(const ReferenceHolder& rhs) MW_NOEXCEPT :
                pImpl(rhs.pImpl) {}

            ReferenceHolder& operator= (const ReferenceHolder& rhs) MW_NOEXCEPT {                            
                pImpl = rhs.pImpl;
                return *this;
            }

        };
    }
}

#endif
