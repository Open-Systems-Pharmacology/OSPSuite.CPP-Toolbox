/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_EXCEPTION_IMPL_HPP
#define ENGINE_EXCEPTION_IMPL_HPP

#include <vector>
#include <streambuf>
#include <memory>
#include <future>
#include "../engine_exception.hpp"
#include <string>
#include "../engine_util.hpp"

#if defined(_WIN32 ) 
#define NOEXCEPT throw() 
#else
#define NOEXCEPT noexcept
#endif


namespace matlab {
    namespace engine { 

        inline Exception::Exception() :std::exception() {}

        inline Exception::Exception(const std::string& msg) : message(msg) {}

        inline Exception::~Exception() {}

        inline Exception& Exception::operator=(const Exception& rhs) {
            message = rhs.message;
            return *this;
        }

        inline const char* Exception::what() const NOEXCEPT { return message.c_str(); }

        inline EngineException::EngineException() {}

        inline EngineException::EngineException(const std::string& msg) : Exception(msg) {}

        inline MATLABNotAvailableException::MATLABNotAvailableException() : Exception() {}

        inline MATLABNotAvailableException::MATLABNotAvailableException(const std::string& msg) : Exception(msg) {}

        inline CancelledException::CancelledException() : Exception() {}

        inline CancelledException::CancelledException(const std::string& msg) : Exception(msg) {}

        inline InterruptedException::InterruptedException() : Exception() {}

        inline InterruptedException::InterruptedException(const std::string& msg) : Exception(msg) {}

        inline MATLABException::MATLABException() : Exception() {}

        inline MATLABException::MATLABException(const std::string& msg) : Exception(msg) {}

        inline MATLABException::MATLABException(const std::string& id, const std::basic_string<char16_t>& txt):Exception(convertUTF16StringToUTF8String(txt)), errorID(id), errorText(txt) {}

        inline MATLABException::MATLABException(const MATLABException& rhs) : Exception(convertUTF16StringToUTF8String(rhs.errorText)), errorID(rhs.errorID), errorText(rhs.errorText) {}

        inline MATLABException& MATLABException::operator=(const MATLABException& rhs) {
            Exception::operator=(convertUTF16StringToUTF8String(rhs.errorText));
            errorID = rhs.errorID;
            errorText = rhs.errorText;
            return *this;
        }

        inline std::string MATLABException::getMessageID() const { return errorID; }

        inline std::basic_string<char16_t> MATLABException::getMessageText() const { return errorText; }

        inline MATLABSyntaxException::MATLABSyntaxException(const std::string& msg) : MATLABException(msg) {}

        inline MATLABSyntaxException::MATLABSyntaxException(const std::string& id, const std::basic_string<char16_t>& txt) : MATLABException(id, txt) {}

        inline MATLABExecutionException::MATLABExecutionException() : MATLABException() {}

        inline MATLABExecutionException::MATLABExecutionException(const std::string& msg, const std::vector<StackFrame>& trace) : MATLABException(msg), stackTrace(trace) {}

        inline MATLABExecutionException::MATLABExecutionException(const std::string& id, const std::basic_string<char16_t>& txt, const std::vector<StackFrame>& trace, const std::vector<MATLABExecutionException>& cause) : MATLABException(id, txt), stackTrace(trace), errorCause(cause) {
        }

        inline std::vector<StackFrame> MATLABExecutionException::getStackTrace() const { return stackTrace; }

        /**
        * Get the cause of the MATLAB exception
        */
        inline std::vector<MATLABExecutionException> MATLABExecutionException::getCause() const { return errorCause; }

        /**
        * Set the cause of the MATLAB exception
        */
        inline void MATLABExecutionException::setCause(const std::vector<MATLABExecutionException>& cause) { errorCause = cause; }

        inline TypeConversionException::TypeConversionException() : Exception() {}

        inline TypeConversionException::TypeConversionException(const std::string& msg) : Exception(msg) {}

        inline StackFrame::StackFrame() {}

        inline StackFrame::StackFrame(const std::basic_string<char16_t>& file, const std::basic_string<char16_t>& func, uint64_t line) :fileName(file), funcName(func), lineNumber(line) {}

        inline std::basic_string<char16_t> StackFrame::getFileName() const {
            return fileName;
        }

        inline std::basic_string<char16_t> StackFrame::getFunctionName() const {
            return funcName;
        }

        inline uint64_t StackFrame::getLineNumber() const {
            return lineNumber;
        }
    }
}


#endif /* ENGINE_EXCEPTION_IMPL_HPP */