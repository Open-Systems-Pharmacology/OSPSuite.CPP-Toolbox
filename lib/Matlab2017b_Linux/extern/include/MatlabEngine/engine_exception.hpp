/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_EXCEPTION_HPP
#define ENGINE_EXCEPTION_HPP

#include <vector>
#include <streambuf>
#include <memory>
#include <future>
#include <exception>
#include <MatlabDataArray/Exception.hpp>

#if defined(_WIN32 ) 
#define NOEXCEPT throw() 
#else
#define NOEXCEPT noexcept
#endif

namespace matlab {

    namespace engine {

        class Exception : public std::exception, public matlab::Exception {
        public:
            Exception();
            
            /**
            * Constructor that accepts an error message
            */
            Exception(const std::string& msg);

            virtual ~Exception();

            Exception& operator=(const Exception& rhs);

            /**
            * Returns the error message
            */
            virtual const char* what() const NOEXCEPT;
        private:
            std::string message;
        };

        class EngineException final: public Exception {
        public:
            EngineException();
            
            /**
            * Constructor that accepts an error message
            */
            EngineException(const std::string& msg);
        };

        class MATLABNotAvailableException final : public Exception {
        public:
            MATLABNotAvailableException();

            /**
            * Constructor that accepts an error message
            */
            MATLABNotAvailableException(const std::string& msg);
        };

        class CancelledException final: public Exception {
        public:
            CancelledException();

            /**
            * Constructor that accepts an error message
            */
            CancelledException(const std::string& msg);
        };

        class InterruptedException final: public Exception {
        public:
            InterruptedException();

            /**
            * Constructor that accepts an error message
            */
            InterruptedException(const std::string& msg);
        };

        class StackFrame final {
        public:
            StackFrame();
            
            /**
            * Constructor that accepts a file name, function name and line number
            */
            StackFrame(const std::basic_string<char16_t>& file, const std::basic_string<char16_t>& func, uint64_t line);

            /**
            * Get the file name
            */
            std::basic_string<char16_t> getFileName() const;

            /**
            * Get the function name
            */
            std::basic_string<char16_t> getFunctionName() const;

            /**
            * Get the line number
            */
            uint64_t getLineNumber() const;
        private:
            std::basic_string<char16_t> fileName;
            std::basic_string<char16_t> funcName;
            uint64_t lineNumber;
        };

        class MATLABException : public Exception {
        public:
            MATLABException();

            /**
            * Constructor that accepts an error message
            */
            MATLABException(const std::string& msg);

            /**
            * Constructor that accepts an error message ID and an error message
            */
            MATLABException(const std::string& id, const std::basic_string<char16_t>& txt);

            /**
            * Get the message ID
            */
            std::string getMessageID() const;

            /**
            * Get the error message
            */
            std::basic_string<char16_t> getMessageText() const;

            virtual ~MATLABException() {};

            /**
            * Copy constructor
            */
            MATLABException(const MATLABException& rhs);

            /**
            * Assignment operator
            */
            MATLABException& operator=(const MATLABException& rhs);
        private:
            std::string errorID;
            std::basic_string<char16_t> errorText;
        };

        class MATLABSyntaxException final : public MATLABException {
        public:
            /**
            * Constructor that accepts an error message
            */
            MATLABSyntaxException(const std::string& msg);

            /**
            * Constructor that accepts an error message ID and an error message
            */
            MATLABSyntaxException(const std::string& id, const std::basic_string<char16_t>& txt);
        };

        class MATLABExecutionException final : public MATLABException {
        public:
            MATLABExecutionException();

            /**
            * Constructor that accepts an error message ID and an error message
            */
            MATLABExecutionException(const std::string& msg, const std::vector<StackFrame>& trace);

            /**
            * Constructor that accepts an error message ID and an error message
            */
            MATLABExecutionException(const std::string& id, const std::basic_string<char16_t>& txt, const std::vector<StackFrame>& trace, const std::vector<MATLABExecutionException>& cause);

            /**
            * Get the stack trace
            */
            std::vector<StackFrame> getStackTrace() const;

            /**
            * Get the cause of the MATLAB exception
            */
            std::vector<MATLABExecutionException> getCause() const;

            /**
            * Set the cause of the MATLAB exception
            */
            void setCause(const std::vector<MATLABExecutionException>& cause);
        private:
            std::vector<StackFrame> stackTrace;
            std::vector<MATLABExecutionException> errorCause;
        };

        class TypeConversionException final : public Exception {
        public:
            TypeConversionException();

            /**
            * Constructor that accepts an error message
            */
            TypeConversionException(const std::string& msg);
        };
    }
}

#endif /* ENGINE_EXCEPTION_HPP */