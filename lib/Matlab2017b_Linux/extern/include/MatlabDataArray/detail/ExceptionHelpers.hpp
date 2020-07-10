/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_EXTDATA_EXCEPTION_HELPER_FUNCTIONS_HPP
#define MATLAB_EXTDATA_EXCEPTION_HELPER_FUNCTIONS_HPP

#include "ExceptionType.hpp"

namespace matlab {
    namespace data {
        namespace detail {
            inline void throwIfError(int type) {
                switch (ExceptionType(type)) {
                  case ExceptionType::NoException:
                    break;
                  case ExceptionType::InvalidFieldName:
                    throw InvalidFieldNameException("Field does not exist in this struct");
                  case ExceptionType::InvalidDataType:
                    throw TypeMismatchException("Data type mismatch");
                  case ExceptionType::InvalidArrayType:
                    throw InvalidArrayTypeException("Array type mismatch");
                  case ExceptionType::MustSpecifyClassName:
                    throw MustSpecifyClassNameException("Must specify a class name");
                  case ExceptionType::WrongNumberOfEnumsSupplied:
                    throw WrongNumberOfEnumsSuppliedException("Wrong number of enums supplied");
                  case ExceptionType::OutOfMemory:
                    throw matlab::OutOfMemoryException("Not enough memory available to support the request");
                  case ExceptionType::TooManyIndicesProvided:
                    throw TooManyIndicesProvidedException("Too many indices provided");
                  case ExceptionType::NotEnoughIndicesProvided:
                    throw NotEnoughIndicesProvidedException("Not enough indices provided");
                  case ExceptionType::StringIndexMustBeLast:
                    throw StringIndexMustBeLastException("String index must be last");
                  case ExceptionType::StringIndexNotValid:
                    throw StringIndexNotValidException("String index not valid for this array");
                  case ExceptionType::CanOnlyUseOneStringIndex:
                    throw CanOnlyUseOneStringIndexException("Can only provide one string index");
                  case ExceptionType::InvalidArrayIndex:
                    throw InvalidArrayIndexException("Array index is invalid");
                  case ExceptionType::CantAssignArrayToThisArray:
                    throw CantAssignArrayToThisArrayException("Can't assign an Array to this array");
                  case ExceptionType::DuplicateFieldNameInStructArray:
                    throw DuplicateFieldNameInStructArrayException("Duplicate field names not allowed");
                  case ExceptionType::NonAsciiCharInRequestedAsciiOutput:
                    throw NonAsciiCharInRequestedAsciiOutputException("Can't create a std::string from this data");
                  default:
                    throw SystemErrorException("Error?");
                }
            }
        }
    }
}

#endif
