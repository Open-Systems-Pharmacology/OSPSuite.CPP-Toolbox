/* Copyright 2017 The MathWorks, Inc. */

#ifndef TASK_REFERENCE_IMPL_HPP
#define TASK_REFERENCE_IMPL_HPP

#include <locale>
#include <string>
#ifdef _MSC_VER
#include <codecvt>
#endif
#include "../cpp_engine_api.hpp"
#include "../task_reference.hpp"

namespace matlab {

    namespace engine {

        inline TaskReference::TaskReference() {}
        inline TaskReference::TaskReference(std::function<bool(uintptr_t, bool)>&& cancel):handle(0), cancelImpl(std::move(cancel)){};
        inline TaskReference::TaskReference(uintptr_t a_handle, std::function<bool(uintptr_t, bool)>&& cancel) :handle(a_handle), cancelImpl(std::move(cancel)) {};
        inline uintptr_t TaskReference::getHandle() const { return handle;  }
        inline TaskReference::~TaskReference(){
            if(handle) engine_destroy_task_handle(handle);
        }
        inline bool TaskReference::cancel(bool allowInterrupt) {
           return cancelImpl(handle, allowInterrupt);
        }
    }
}

#endif //TASK_REFERENCE_IMPL_HPP