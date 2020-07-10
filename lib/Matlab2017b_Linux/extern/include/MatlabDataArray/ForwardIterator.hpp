/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_DATA_FORWARD_ITERATOR_HPP
#define MATLAB_DATA_FORWARD_ITERATOR_HPP

#include <memory>
#include <iterator>
#include <functional>
        
#include "detail/forward_iterator_interface.hpp"

namespace matlab {
    namespace data {
        namespace detail {
            class Access;
            class RefCounted;
            class ForwardIteratorImpl;
        }

        template<typename T>
        class ForwardIterator {
          public:

            typedef std::ptrdiff_t difference_type;
            typedef T value_type;
            typedef T* pointer;
            typedef T reference;
            typedef std::forward_iterator_tag iterator_category;

            /**
             * ForwardIterator copy constructor. Guaranteed not to throw.
             *
             * @param rhs ForwardIterator object to be copied.
             * @throw none.
             */
            ForwardIterator(const ForwardIterator<T> &rhs) MW_NOEXCEPT :
                pImpl(std::shared_ptr<detail::ForwardIteratorImpl>
                      (forward_iterator_clone(rhs.pImpl.get()), 
                       [](detail::ForwardIteratorImpl* ptr) {
                          forward_iterator_destroy_impl(ptr);
                      })) {}
            
            /**
             * ForwardIterator operator=.
             *
             * @param rhs ForwardIterator object to assign to this instance.
             * @return the updated instance.
             * @throw none
             */
            ForwardIterator<T>& operator=(const ForwardIterator<T> &rhs) MW_NOEXCEPT {
                pImpl = std::shared_ptr<detail::ForwardIteratorImpl>
                    (forward_iterator_clone(rhs.pImpl.get()), 
                     [](detail::ForwardIteratorImpl* ptr) {
                        forward_iterator_destroy_impl(ptr);
                    });
                fRef = rhs.fRef;
                return *this;
            }
 
            /**
             * Pre-increment 
             *
             * @return ForwardIterator<T>& - reference to this.
             * @throw none
             */
            ForwardIterator<T>& operator++() MW_NOEXCEPT {
                forward_iterator_plus_plus(pImpl.get());
                return *this;
            }

            /**
             * Post-increment
             *
             * @return ForwardIterator<T> - newly created ForwardIterator
             * @throw none
             */
            ForwardIterator<T> operator++(int) MW_NOEXCEPT {
                ForwardIterator<T> temp(*this);
                forward_iterator_plus_plus(pImpl.get());
                return temp;
            }

            /**
             * operator== determine if 2 iterators are pointing to the same element
             * @param rhs - the iterator to be compared with
             * @return bool - true if the iterators are pointing to the same element
             * @throw none
             */
            bool operator==(const ForwardIterator<T> &rhs) const MW_NOEXCEPT {
                return forward_iterator_equal(pImpl.get(), rhs.pImpl.get());
            }

            /**
             * operator!= determine if 2 iterators are pointing to the same element
             * @param rhs - the iterator to be compared with
             * @return bool - true if this iterator points to a different element
             * @throw none
             */
            bool operator!=(const ForwardIterator<T> &rhs) const MW_NOEXCEPT {
                return !forward_iterator_equal(pImpl.get(), rhs.pImpl.get());
            }

            /**
             * operator* return a shared copy of the element that the iterator is pointing to
             *
             * @return T - shared copy of the element that the iterator is pointing to
             * @throw none
             */
            reference operator*() const MW_NOEXCEPT {
                detail::RefCounted* value = nullptr;
                forward_iterator_get_ref(pImpl.get(), &value);
                return detail::Access::createObj<reference>(value);
            }

            /**
             * Get a pointer the element pointed to by this iterator
             *
             * @return pointer - the element
             * @throw none
             */
            pointer operator->() const MW_NOEXCEPT {
                detail::RefCounted* value = nullptr;
                forward_iterator_get_ref(pImpl.get(), &value);
                fRef = detail::Access::createObj<reference>(value);
                return &fRef;
            }
            
            /**
             * Get a reference using a linear index
             *
             * @return reference - the element
             * @throw none
             */
            reference operator[](const size_t& rhs) const MW_NOEXCEPT {
                ForwardIterator<T> temp(*this);
                forward_iterator_increment(temp.pImpl.get(), rhs);
                return *temp;
            }

          private:
            friend class detail::Access;
            
            std::shared_ptr<detail::ForwardIteratorImpl> pImpl;   
            mutable std::remove_const<reference> fRef;
           
            ForwardIterator(detail::ForwardIteratorImpl *impl) :
                pImpl(std::shared_ptr<detail::ForwardIteratorImpl>
                      (impl, 
                       [](detail::ForwardIteratorImpl* ptr) {
                          forward_iterator_destroy_impl(ptr);
                      })) {}
        };
    }
}


#endif
