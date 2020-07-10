/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_UTIL_HPP
#define ENGINE_UTIL_HPP

#include <locale>
#include <string>
#ifdef _MSC_VER
#include <codecvt>
#endif

#include "cpp_engine_api.hpp"

namespace detail {
    enum class EngineExceptionType {
        CANCELLED = 0,
        INTERRUPTED = 1,
        EXECUTION = 2,
        SYNTAX = 3,
        OTHER = 4,
        STOPPED = 5
    };

    enum class ConvertUTFStringExceptionType {
        CONVERT_ERROR = 1,
        OUT_OF_MEMORY = 2
    };
}


namespace matlab {

    namespace engine {

        typedef std::basic_string<char16_t> String;
        typedef std::basic_streambuf<char16_t> StreamBuffer;
        class MATLABEngine;
        template <class T>
        class FutureResult;

        enum class WorkspaceType {
            BASE = 0,
            GLOBAL = 1
        };

        std::basic_string<char16_t> convertUTF8StringToUTF16String(const std::string &utf8string);

        std::string convertUTF16StringToUTF8String(const std::basic_string<char16_t> &utf16string);
    }
}

#endif //ENGINE_UTIL_HPP