/* Copyright 2017 The MathWorks, Inc. */

#ifndef MATLAB_ENGINE_IMPL_HPP
#define MATLAB_ENGINE_IMPL_HPP

#include <vector>
#include <streambuf>
#include <memory>
#include <future>
#include <complex>

#include <MatlabDataArray/detail/HelperFunctions.hpp>
#include "MatlabDataArray/StructArray.hpp"
#include "MatlabDataArray/StructRef.hpp"
#include "MatlabDataArray/Struct.hpp"
#include "MatlabDataArray/CharArray.hpp"
#include "MatlabDataArray/TypedArray.hpp"
#include "MatlabDataArray/Reference.hpp"
#include "../matlab_engine.hpp"
#include "../engine_exception.hpp"

namespace {

    inline matlab::engine::MATLABExecutionException createMATLABExecutionException(const matlab::data::StructArray& mException);

    inline std::vector<matlab::engine::MATLABExecutionException> createCause(const matlab::data::CellArray& cause) {
        size_t nCauses = cause.getNumberOfElements();
        std::vector<matlab::engine::MATLABExecutionException> causes(nCauses);
        for (size_t i = 0; i < nCauses; i++) {
            matlab::data::Array exRef = cause[i];
            matlab::data::StructArray ex(exRef);
            causes[i] = createMATLABExecutionException(ex);
        }
        return causes;
    }

    inline std::vector<matlab::engine::StackFrame> createStackTrace(const matlab::data::StructArray& stack) {
        size_t nFrames = stack.getNumberOfElements();
        std::vector<matlab::engine::StackFrame> stackFrames(nFrames);

        for (size_t i = 0; i < nFrames; i++) {
            matlab::data::Array fileRef = stack[i]["File"];
            matlab::data::CharArray fileStr(fileRef);

            matlab::data::Array nameRef = stack[i]["Name"];
            matlab::data::CharArray nameStr(nameRef);

            matlab::data::Array lineRef = stack[i]["Line"];
            double line = lineRef[0];
            stackFrames[i] = matlab::engine::StackFrame(fileStr.toUTF16(), nameStr.toUTF16(), uint32_t(line));
        }

        return stackFrames;
    }

    inline matlab::engine::MATLABExecutionException createMATLABExecutionException(const matlab::data::StructArray& mException) {
        matlab::data::Array idRef = mException[0][std::string("identifier")];
        matlab::data::CharArray id(idRef);
        matlab::data::Array messageRef = mException[0][std::string("message")];
        matlab::data::CharArray message(messageRef);

        matlab::data::Array stackRef = mException[0][std::string("stack")];
        matlab::data::StructArray stack(stackRef);

        matlab::data::Array causeRef = mException[0][std::string("cause")];
        matlab::data::CellArray cause(causeRef);

        std::vector<matlab::engine::MATLABExecutionException> meCause = createCause(cause);
        std::vector<matlab::engine::StackFrame> meStack = createStackTrace(stack);
        return matlab::engine::MATLABExecutionException(id.toAscii(), message.toUTF16(), meStack, meCause);
    }

    inline matlab::engine::MATLABSyntaxException createMATLABSyntaxException(const matlab::data::StructArray& mException) {
        matlab::data::Array idRef = mException[0][std::string("identifier")];
        matlab::data::CharArray id(idRef);
        matlab::data::Array messageRef = mException[0][std::string("message")];
        matlab::data::CharArray message(messageRef);

        return matlab::engine::MATLABSyntaxException(id.toAscii(), message.toUTF16());
    }

    template<typename T>
    inline void set_promise_exception(void *p, size_t excTypeNumber, const void* msg) {
        std::promise<T>* prom = reinterpret_cast<std::promise<T>*>(p);
        ::detail::EngineExceptionType excType = static_cast<::detail::EngineExceptionType>(excTypeNumber);
        switch (excType) {
        case ::detail::EngineExceptionType::CANCELLED:{
            const char* message = reinterpret_cast<const char*>(msg);
            matlab::engine::CancelledException exception(message);
            prom->set_exception(std::make_exception_ptr<matlab::engine::CancelledException>(exception));
            break;
        }
        case ::detail::EngineExceptionType::INTERRUPTED: {
            const char* message = reinterpret_cast<const char*>(msg);
            matlab::engine::InterruptedException exception(message);
            prom->set_exception(std::make_exception_ptr<matlab::engine::InterruptedException>(exception));
            break;
        }
        case ::detail::EngineExceptionType::EXECUTION:
        case ::detail::EngineExceptionType::SYNTAX: {
            matlab::data::impl::ArrayImpl* exceptionImpl = const_cast<matlab::data::impl::ArrayImpl*>(reinterpret_cast<const matlab::data::impl::ArrayImpl*>(msg));
            matlab::data::Array mdaException = matlab::data::detail::Access::createObj<matlab::data::Array>(exceptionImpl);
            matlab::data::StructArray mException(mdaException);
            if (excType == ::detail::EngineExceptionType::SYNTAX) {
                matlab::engine::MATLABSyntaxException exception = createMATLABSyntaxException(mException);
                prom->set_exception(std::make_exception_ptr<matlab::engine::MATLABSyntaxException>(exception));
            }
            else {
                matlab::engine::MATLABExecutionException exception(createMATLABExecutionException(mException));
                prom->set_exception(std::make_exception_ptr<matlab::engine::MATLABExecutionException>(exception));
            }
            break;
        }
        case ::detail::EngineExceptionType::OTHER: {
            const char* message = reinterpret_cast<const char*>(msg);
            matlab::engine::EngineException exception(message);
            prom->set_exception(std::make_exception_ptr<matlab::engine::EngineException>(exception));
            break;
        }
        case ::detail::EngineExceptionType::STOPPED: {
            const char* message = reinterpret_cast<const char*>(msg);
            matlab::engine::MATLABNotAvailableException exception(message);
            prom->set_exception(std::make_exception_ptr<matlab::engine::MATLABNotAvailableException>(exception));
            break;
        }
        break;
        }
        delete prom;
    }
}

namespace matlab {
    namespace engine {

        inline MATLABEngine::MATLABEngine(uint64_t handle) : matlabHandle(handle) {
        }

        inline void set_eval_promise_data(void *p) {
            std::promise<void>* prom = reinterpret_cast<std::promise<void>*>(p);
            prom->set_value();
            delete prom;
        }
        
        inline void set_eval_promise_exception(void *p, size_t excTypeNumber, const void* msg) {
            set_promise_exception<void>(p, excTypeNumber, msg);
        }

        inline void set_feval_promise_data(void *p, size_t nlhs, bool straight, matlab::data::impl::ArrayImpl** plhs) {

            if (nlhs == 0 && straight) {
                std::promise<void>* prom = reinterpret_cast<std::promise<void>*>(p);
                prom->set_value();
                delete prom;
                return;
            }

            if (nlhs == 1 && straight) {
                std::promise<matlab::data::Array>* prom = reinterpret_cast<std::promise<matlab::data::Array>*>(p);
                matlab::data::Array v_ = matlab::data::detail::Access::createObj<matlab::data::Array>(plhs[0]);
                prom->set_value(v_);
                delete prom;
                return;
            }

            std::promise<std::vector<matlab::data::Array> >* prom = reinterpret_cast<std::promise<std::vector<matlab::data::Array> >*>(p);
            std::vector<matlab::data::Array> result;
            for (size_t i = 0; i < nlhs; i++) {
                matlab::data::Array v_ = matlab::data::detail::Access::createObj<matlab::data::Array>(plhs[i]);
                result.push_back(v_);
            }
            prom->set_value(result);
            delete prom;
        }

        template<class T>
        void set_exception(T p, std::exception_ptr e) {
            p->set_exception(e);
        }

        inline void set_feval_promise_exception(void *p, size_t nlhs, bool straight, size_t excTypeNumber, const void* msg) {
            if (nlhs == 0 && straight) {
                set_promise_exception<void>(p, excTypeNumber, msg);
            }
            else if (nlhs == 1 && straight) {
                set_promise_exception<matlab::data::Array>(p, excTypeNumber, msg);
            }
            else {
                set_promise_exception<std::vector<matlab::data::Array>>(p, excTypeNumber, msg);
            }
        }
        
        inline std::vector<matlab::data::Array> MATLABEngine::feval(const String &function,
                                                             const size_t nlhs,
                                                             const std::vector<matlab::data::Array> &args,
                                                             const std::shared_ptr<StreamBuffer>& output,
                                                             const std::shared_ptr<StreamBuffer>& error) {
            return fevalAsync(function, nlhs, args, output, error).get();
        }

        
        inline matlab::data::Array MATLABEngine::feval(const String &function,
            const std::vector<matlab::data::Array> &args,
            const std::shared_ptr<StreamBuffer>& output,
            const std::shared_ptr<StreamBuffer>& error) {
            FutureResult<matlab::data::Array> future = fevalAsync(function, args, output, error);
            return future.get();
        }

        inline matlab::data::Array MATLABEngine::feval(const String &function,
        const matlab::data::Array &arg,
        const std::shared_ptr<StreamBuffer>& output,
        const std::shared_ptr<StreamBuffer>& error) {
            FutureResult<matlab::data::Array> future = fevalAsync(function, arg, output, error);
            return future.get();
        }

        template<class ReturnType, typename...RhsArgs>
        ReturnType MATLABEngine::feval(const String &function,
        const std::shared_ptr<StreamBuffer>& output,
        const std::shared_ptr<StreamBuffer>& error,
        RhsArgs&&... rhsArgs
        ) {
            return fevalAsync<ReturnType>(function, output, error, std::forward<RhsArgs>(rhsArgs)...).get();
        }
        

        namespace detail {
            template<typename T>
            inline void validateTIsSupported() {
                using U = typename std::remove_cv<typename std::remove_reference<T>::type>::type;
                static_assert(
                    std::is_same<U, bool>::value
                    || std::is_same<U, int>::value
                    || std::is_same<U, int8_t>::value
                    || std::is_same<U, int16_t>::value
                    || std::is_same<U, int32_t>::value
                    || std::is_same<U, int64_t>::value
                    || std::is_same<U, uint8_t>::value
                    || std::is_same<U, uint16_t>::value
                    || std::is_same<U, uint32_t>::value
                    || std::is_same<U, uint64_t>::value
                    || std::is_same<U, float>::value
                    || std::is_same<U, double>::value, "Attempted to use unsupported types.");
            }

            template<class T>
            matlab::data::Array createRhs(matlab::data::ArrayFactory& factory, T&& value) {
                validateTIsSupported<T>();
                using U = typename std::remove_cv<typename std::remove_reference<T>::type>::type;
                return factory.createArray<U>({ 1, 1 }, {value});
            }

            template<typename T, typename A>
            matlab::data::Array createRhs(matlab::data::ArrayFactory& factory, std::vector <T, A>&& value) {
                validateTIsSupported<T>();
                return factory.createArray({ 1, value.size() }, value.begin(), value.end());
            }

            template <std::size_t ...Ints>
            struct index_sequence {
                using value_type = std::size_t;
                static std::size_t size() { return sizeof...(Ints); }
            };

            template<std::size_t N, std::size_t... Values>
            struct make_index_sequence_impl {
                using type = typename make_index_sequence_impl<N - 1, Values..., sizeof...(Values)>::type;
            };

            template<std::size_t... Values>
            struct make_index_sequence_impl < 0, Values... > {
                using type = index_sequence < Values... > ;
            };

            template<std::size_t N>
            using make_index_sequence = typename make_index_sequence_impl<N>::type;

            template<typename T>
            struct createLhs {
                static const size_t nlhs = 1;

                T operator()(std::vector<matlab::data::Array>&& lhs) const {
                    if (lhs.empty()) {
                        throw matlab::engine::TypeConversionException("The result is empty.");
                    }

                    T value;
                    try {
                        value = (*this)(matlab::data::TypedArray<T>(std::move(lhs.front())));
                    }
                    catch (const std::exception& e) {
                        throw matlab::engine::TypeConversionException(e.what());
                    }
                    return value;
                }

                T operator()(matlab::data::TypedArray<T> lhs) const {
                    validateTIsSupported<T>();
                    auto const begin = lhs.begin();
                    auto const end = lhs.end();
                    if (begin == end) {
                        throw matlab::engine::TypeConversionException("The result is empty.");
                    }
                    return *begin;
                }
            };

            template<>
            struct createLhs < void > {
                static const size_t nlhs = 0;
                void operator()(std::vector<matlab::data::Array>&& lhs) const {}
            };

            template<typename... TupleTypes>
            struct createLhs < std::tuple<TupleTypes...> > {
                static const size_t nlhs = sizeof...(TupleTypes);
                using T = std::tuple < TupleTypes... > ;

                T operator()(std::vector<matlab::data::Array>&& lhs) const {
                    //we are not validating the LHS here as it can be any combinations of types for std::tuple.
                    if (lhs.size() < sizeof...(TupleTypes)) { throw std::runtime_error(""); }
                    return (*this)(std::move(lhs), detail::make_index_sequence<sizeof...(TupleTypes)>());
                }

            private:
                template<size_t Index>
                using TupleElement = typename std::remove_cv<typename std::remove_reference<typename std::tuple_element<Index, std::tuple<TupleTypes...> >::type>::type>::type;

                template<size_t... IndexList>
                std::tuple <TupleTypes...> operator()(std::vector<matlab::data::Array>&& lhs, detail::index_sequence<IndexList...>) const {
                    return std::tuple <TupleTypes...>(createLhs<TupleElement<IndexList>>()(std::move(lhs[IndexList]))...);
                }
            };

        }

        template<class ReturnType, typename...RhsArgs>
        ReturnType MATLABEngine::feval(const String &function,
                              RhsArgs&&...rhsArgs
        ) {
            const std::shared_ptr<StreamBuffer> defaultStream;
            auto future = fevalAsync<ReturnType>(function, defaultStream, defaultStream, std::forward<RhsArgs>(rhsArgs)...);
            return future.get();
        }
        
        inline void MATLABEngine::eval(const String &str,
        const std::shared_ptr<StreamBuffer> &output,
        const std::shared_ptr<StreamBuffer> &error
        ) {
            return evalAsync(str, output, error).get();
        }

        inline matlab::data::Array MATLABEngine::getVariable(const String &varName, WorkspaceType workspaceType) {
            return getVariableAsync(varName).get();
        }

        inline void MATLABEngine::setVariable(const String &varName, const matlab::data::Array& var, WorkspaceType workspaceType) {
            return setVariableAsync(varName, var, workspaceType).get();
        }

        inline matlab::data::Array MATLABEngine::getProperty(const matlab::data::Array &object, const String &propertyName) {
            return getPropertyAsync(object, propertyName).get();
        }

        inline void MATLABEngine::setProperty(matlab::data::Array &object, const String &propertyName, const matlab::data::Array &property) {
            return setPropertyAsync(object, propertyName, property).get();
        }

        inline MATLABEngine::~MATLABEngine() {
            if (matlabHandle != 0) {
                engine_terminate_out_of_process_matlab(matlabHandle);
                matlabHandle = 0;
            }
        }

        inline FutureResult<std::vector<matlab::data::Array> > MATLABEngine::fevalAsync(const String &function,
        const size_t nlhs,
        const std::vector<matlab::data::Array> &args,
        const std::shared_ptr<StreamBuffer> &output,
        const std::shared_ptr<StreamBuffer> &error
        ) {
            size_t nrhs = args.size();
            matlab::data::impl::ArrayImpl** argsImpl = new matlab::data::impl::ArrayImpl*[nrhs];
            size_t i = 0;
            for (auto e : args) {
                argsImpl[i++] = matlab::data::detail::Access::getImpl<matlab::data::impl::ArrayImpl>(e);
            }
            std::promise<std::vector<matlab::data::Array> >* p = new std::promise<std::vector<matlab::data::Array> >();
            std::future<std::vector<matlab::data::Array> > f = p->get_future();

            void* output_ = output ? new std::shared_ptr<StreamBuffer>(std::move(output)) : nullptr;
            void* error_ = error ? new std::shared_ptr<StreamBuffer>(std::move(error)) : nullptr;

            uintptr_t handle = engine_feval_with_completion(matlabHandle, convertUTF16StringToUTF8String(function).c_str(), nlhs, false, argsImpl, nrhs, &set_feval_promise_data, &set_feval_promise_exception, p, output_, error_, &writeToStreamBuffer, &deleteStreamBufferImpl);

            return FutureResult<std::vector<matlab::data::Array>>(std::move(f), std::make_shared<TaskReference>(handle, &engine_cancel_feval_with_completion));
        }

        inline FutureResult<matlab::data::Array> MATLABEngine::fevalAsync(const String &function,
        const std::vector<matlab::data::Array> &args,
        const std::shared_ptr<StreamBuffer> &output,
        const std::shared_ptr<StreamBuffer> &error
        ) {
            size_t nrhs = args.size();
            matlab::data::impl::ArrayImpl** argsImpl = new matlab::data::impl::ArrayImpl*[nrhs];
            size_t i = 0;
            for (auto e : args) {
                argsImpl[i++] = matlab::data::detail::Access::getImpl<matlab::data::impl::ArrayImpl>(e);
            }
            std::promise<matlab::data::Array>* p = new std::promise<matlab::data::Array>();
            std::future<matlab::data::Array> f = p->get_future();

            void* output_ = output ? new std::shared_ptr<StreamBuffer>(std::move(output)) : nullptr;
            void* error_ = error ? new std::shared_ptr<StreamBuffer>(std::move(error)) : nullptr;

            uintptr_t handle = engine_feval_with_completion(matlabHandle, convertUTF16StringToUTF8String(function).c_str(), 1, true, argsImpl, nrhs, &set_feval_promise_data, &set_feval_promise_exception, p, output_, error_, &writeToStreamBuffer, &deleteStreamBufferImpl);

            return FutureResult<matlab::data::Array>(std::move(f), std::make_shared<TaskReference>(handle, &engine_cancel_feval_with_completion));
        }

        inline FutureResult<matlab::data::Array> MATLABEngine::fevalAsync(const String &function,
        const matlab::data::Array &arg,
        const std::shared_ptr<StreamBuffer> &output,
        const std::shared_ptr<StreamBuffer> &error
        ) {
            return fevalAsync(function, std::vector<matlab::data::Array>({ arg }), output, error);
        }

        template<class ReturnType, typename... RhsArgs>
        FutureResult<ReturnType> MATLABEngine::fevalAsync(const String &function,
            const std::shared_ptr<StreamBuffer> &output,
            const std::shared_ptr<StreamBuffer> &error,
            RhsArgs&&... rhsArgs
            ) {
            matlab::data::ArrayFactory factory;
            std::vector<matlab::data::Array> rhsList({
                detail::createRhs(factory, std::forward<RhsArgs>(rhsArgs))...
            });

            auto const nlhs = detail::createLhs<ReturnType>::nlhs;

            size_t nrhs = rhsList.size();
            matlab::data::impl::ArrayImpl** argsImpl = new matlab::data::impl::ArrayImpl*[nrhs];
            size_t i = 0;
            for (auto e : rhsList) {
                argsImpl[i++] = matlab::data::detail::Access::getImpl<matlab::data::impl::ArrayImpl>(e);
            }

            FutureResult<std::vector<matlab::data::Array>> f = fevalAsync(function, nlhs, rhsList, output, error);

            // c++11 lambdas do not correctly handle move operations...
            // when c++14 is available, this could should be:
            // auto converToResulType = [copyableF = std::move(f)]()->ReturnType { ....... };
            auto copyableF = std::make_shared<FutureResult<std::vector<matlab::data::Array>>>(std::move(f));
            auto convertToResultType = [copyableF]() ->ReturnType {
                std::vector<matlab::data::Array> vec = copyableF->get();
                detail::createLhs<ReturnType> lhsFactory;
                return lhsFactory(std::move(vec));
            };

            std::future<ReturnType> future = std::async(std::launch::deferred, std::move(convertToResultType));
            return FutureResult<ReturnType>(std::move(future), copyableF->getTaskReference());
        }


        template<class ReturnType, typename... RhsArgs>
        FutureResult<ReturnType> MATLABEngine::fevalAsync(const String &function,
            RhsArgs&&... rhsArgs
            ) {
            const std::shared_ptr<StreamBuffer> defaultBuffer;
            return fevalAsync<ReturnType>(function, defaultBuffer, defaultBuffer, std::forward<RhsArgs>(rhsArgs)...);
        }
        
        inline FutureResult<void> MATLABEngine::evalAsync(const String &statement,
            const std::shared_ptr<StreamBuffer> &output,
            const std::shared_ptr<StreamBuffer> &error
            ) {
           
            std::promise<void>* p = new std::promise<void>();
            std::future<void> f = p->get_future();

            uintptr_t* handles;
            std::vector<uintptr_t> handlesV(3);
            
            void* output_ = output? new std::shared_ptr<StreamBuffer>(std::move(output)) : nullptr;
            void* error_ = error? new std::shared_ptr<StreamBuffer>(std::move(error)) : nullptr;

            engine_eval_with_completion(matlabHandle, statement.c_str(), &set_eval_promise_data, &set_eval_promise_exception, p, output_, error_, &writeToStreamBuffer, &deleteStreamBufferImpl, &handles);

            for (size_t i = 0; i < 3; i++) {
                handlesV[i] = handles[i];
            }
            engine_destroy_handles(handles);
            return FutureResult<void>(std::move(f), std::make_shared<TaskReference>(handlesV[0], &engine_cancel_feval_with_completion));
        }

        namespace {
            matlab::engine::String getUString(std::string varName) {
                const char* varName_c = varName.c_str();
                std::basic_string<char16_t> uVarName;
                for (size_t j = 0; j < varName.size(); j++) {
                    uVarName += static_cast<char16_t>(varName_c[j]);
                }
                return uVarName;
            }
        }
        inline FutureResult<matlab::data::Array> MATLABEngine::getVariableAsync(const String &varName, WorkspaceType workspaceType) {
            matlab::data::ArrayFactory factory;
            matlab::data::Array arg = factory.createCharArray(varName);
            return fevalAsync(convertUTF8StringToUTF16String("matlab.internal.engine.getVariable"), arg, nullptr, nullptr);
        }

        inline FutureResult<void> MATLABEngine::setVariableAsync(const String &varName, const matlab::data::Array& var, WorkspaceType workspaceType) {
            matlab::data::ArrayFactory factory;
            std::string base = (workspaceType == WorkspaceType::BASE) ? std::string("base") : std::string("global");
            auto arg1 = factory.createCharArray(base);
            auto arg2 = factory.createCharArray(varName);
            std::vector<matlab::data::Array> args;
            args.push_back(arg2);
            args.push_back(var);
            args.push_back(arg1);
            size_t nrhs = args.size();
            matlab::data::impl::ArrayImpl** argsImpl = new matlab::data::impl::ArrayImpl*[nrhs];
            size_t i = 0;
            for (auto e : args) {
                argsImpl[i++] = matlab::data::detail::Access::getImpl<matlab::data::impl::ArrayImpl>(e);
            }
            std::promise<void>* p = new std::promise<void>();
            std::future<void> f = p->get_future();
            uintptr_t handle = engine_feval_with_completion(matlabHandle, "matlab.internal.engine.setVariable", 0, true, argsImpl, nrhs, &set_feval_promise_data, &set_feval_promise_exception, p, nullptr, nullptr, &writeToStreamBuffer, &deleteStreamBufferImpl);
            return FutureResult<void>(std::move(f), std::make_shared<TaskReference>(handle, &engine_cancel_feval_with_completion));
        }

        inline FutureResult<matlab::data::Array> MATLABEngine::getPropertyAsync(const matlab::data::Array &object, const String &propertyName) {
            matlab::data::ArrayFactory factory;
            auto arg = factory.createCharArray(propertyName);
            std::vector<matlab::data::Array> args;
            args.push_back(object);
            args.push_back(arg);
            return fevalAsync(convertUTF8StringToUTF16String("matlab.internal.engine.getProperty"), args, nullptr, nullptr);
        }

        inline FutureResult<void> MATLABEngine::setPropertyAsync(matlab::data::Array &object, const String &propertyName, const matlab::data::Array &property) {
            matlab::data::ArrayFactory factory;
            auto arg = factory.createCharArray(propertyName);
            std::vector<matlab::data::Array> args;
            args.push_back(object);
            args.push_back(arg);
            args.push_back(property);
            size_t nrhs = args.size();
            FutureResult<void> ret;
            if (object.getType() == matlab::data::ArrayType::VALUE_OBJECT) {
                FutureResult<matlab::data::Array> f = fevalAsync(convertUTF8StringToUTF16String("matlab.internal.engine.setProperty"), args);
                auto copyableF = std::make_shared<FutureResult<matlab::data::Array>>(std::move(f));
                auto convertToResultType = [&object, copyableF]() {
                    matlab::data::Array vec = copyableF->get();
                    object = std::move(vec);
                };
                std::future<void> future = std::async(std::launch::deferred, std::move(convertToResultType));
                ret = FutureResult<void>(std::move(future), copyableF->getTaskReference());
            }
            else if (object.getType() == matlab::data::ArrayType::HANDLE_OBJECT_REF) {
                matlab::data::impl::ArrayImpl** argsImpl = new matlab::data::impl::ArrayImpl*[nrhs];
                size_t i = 0;
                for (auto e : args) {
                    argsImpl[i++] = matlab::data::detail::Access::getImpl<matlab::data::impl::ArrayImpl>(e);
                }
                std::promise<void>* p = new std::promise<void>();
                std::future<void> f = p->get_future();
                uintptr_t handle = engine_feval_with_completion(matlabHandle, "matlab.internal.engine.setProperty", 0, true, argsImpl, nrhs, &set_feval_promise_data, &set_feval_promise_exception, p, nullptr, nullptr, &writeToStreamBuffer, &deleteStreamBufferImpl);

                ret = std::move(FutureResult<void>(std::move(f), std::make_shared<TaskReference>(handle, &engine_cancel_feval_with_completion)));
            }
            else {
                throw EngineException("The input variable is not a MATLAB object.");
            }
            return ret;
        }

        inline void writeToStreamBuffer(void* buffer, const char16_t* stream, size_t n) {
            std::shared_ptr<StreamBuffer>* output = reinterpret_cast<std::shared_ptr<StreamBuffer>*>(buffer);
            output->get()->sputn(stream, n);
        }

        inline void deleteStreamBufferImpl(void* impl) {
            delete static_cast<std::shared_ptr<StreamBuffer>*>(impl);
        }
    }
}

#endif /* MATLAB_ENGINE_IMPL_HPP */
