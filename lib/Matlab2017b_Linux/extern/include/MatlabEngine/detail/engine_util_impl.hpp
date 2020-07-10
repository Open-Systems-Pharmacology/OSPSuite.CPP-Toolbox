/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_UTIL_IMPL_HPP
#define ENGINE_UTIL_IMPL_HPP

#include <locale>
#include <string>
#ifdef _MSC_VER
#include <codecvt>
#endif
#include "../engine_exception.hpp"
#include "../engine_util.hpp"

namespace matlab {

    namespace engine {

        inline std::basic_string<char16_t> convertUTF8StringToUTF16String(const std::string &utf8string) {
            char16_t* name = nullptr;
            std::unique_ptr<char16_t, void(*)(char16_t*)> guard(name, &engine_destroy_utf16);
            size_t errType = 0;
            engine_utf8_to_utf16(utf8string.c_str(), &name, &errType);
            if (errType > 0) {
                ::detail::ConvertUTFStringExceptionType excType = static_cast<::detail::ConvertUTFStringExceptionType>(errType);
                //it could be either out of memory
                if (excType == ::detail::ConvertUTFStringExceptionType::OUT_OF_MEMORY) throw(matlab::OutOfMemoryException("Not enough memory available to support the request."));
                // or std::range_error during the conversion
                if (excType ==  ::detail::ConvertUTFStringExceptionType::CONVERT_ERROR) throw(TypeConversionException("The input cannot be converted to the requested type."));
            }
            String uName(name);
            return uName;
        }

        inline std::string convertUTF16StringToUTF8String(const std::basic_string<char16_t> &utf16string) {
            char* name = nullptr;
            std::unique_ptr<char, void(*)(char*)> guard(name, &engine_destroy_utf8);
            size_t errType = 0;
            engine_utf16_to_utf8(utf16string.c_str(), &name, &errType);
            if (errType > 0) {
                ::detail::ConvertUTFStringExceptionType excType = static_cast<::detail::ConvertUTFStringExceptionType>(errType);
                //it could be either out of memory
                if (excType == ::detail::ConvertUTFStringExceptionType::OUT_OF_MEMORY) throw(matlab::OutOfMemoryException("Not enough memory available to support the request."));
                // or std::range_error during the conversion
                if (excType == ::detail::ConvertUTFStringExceptionType::CONVERT_ERROR) throw(TypeConversionException("The input cannot be converted to the requested type."));
            }
            std::string sName(name);
            return sName;
        }
    }
}

#endif //ENGINE_UTIL_IMPL_HPP