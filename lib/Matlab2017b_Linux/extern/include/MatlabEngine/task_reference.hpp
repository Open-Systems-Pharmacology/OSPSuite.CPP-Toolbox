/* Copyright 2017 The MathWorks, Inc. */

#ifndef TASK_REFERENCE_HPP
#define TASK_REFERENCE_HPP

#include <string>
#include <functional>

namespace matlab {

    namespace engine {

        class TaskReference {
        public:
            TaskReference();

            /**
            * A constructor that accepts a functor for cancel behavior
            * 
            * @param cancel - A functor that used to cancel a task
            */
            TaskReference(std::function<bool(uintptr_t, bool)>&& cancel);

            /**
            * A constructor that accepts a handle and a functor for cancel behavior
            *
            * @param a_handle - The handle of the task
            * @param cancel - A functor that used to cancel a task
            *
            */
            TaskReference(uintptr_t a_handle, std::function<bool(uintptr_t, bool)>&& cancel);

            /**
            * Get the handle of the task
            *
            * @return the handle of the task
            *
            * @throw none
            */
            uintptr_t getHandle() const;

            virtual ~TaskReference();

            /**
            * Cancel the task
            *
            * @param allowInterrupt - Interrupt or not if the task is being processed
            * @return true if the task is cancel-able; false otherwise
            *
            * @throw none
            */
            virtual bool cancel(bool allowInterrupt);
        private:
            TaskReference(TaskReference& handle) = delete;
            TaskReference& operator=(TaskReference&) = delete;
            uintptr_t handle;
            std::function<bool(uintptr_t, bool)> cancelImpl;
            
        };
    }
}



#endif //TASK_REFERENCE_HPP