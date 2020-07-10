/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_HELPER_FUNCTIONS_HPP
#define MATLAB_EXTDATA_HELPER_FUNCTIONS_HPP

#include <memory>

namespace matlab {
    namespace data {
        namespace detail {
            class Access {
            public:
                template<typename ImplT, typename ObjT>
                static ImplT* getImpl(const ObjT& a) {
                    return static_cast<ImplT *>(a.pImpl.get());
                }
                template <typename Obj, typename Impl>
                static Obj createObj(Impl* impl) {
                    return Obj(impl);
                }
                template <typename Obj, typename Impl>
                static Obj createObj(std::shared_ptr<Impl> impl) {
                    return Obj(impl);
                }
            };
        }
    }
}

#endif
