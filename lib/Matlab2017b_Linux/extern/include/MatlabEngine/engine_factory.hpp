/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_FACTORY_HPP
#define ENGINE_FACTORY_HPP

#include <vector>
#include <string>
#include <memory>
#include "engine_util.hpp"

namespace matlab {

    namespace engine {

        class MATLABEngine;

        /**
         * Start MATLAB synchronously with startup options
         *
         * @param options -  Startup options used to launch MATLAB
         * @return std::unique_ptr<MATLABEngine> - A unique pointer of MATLABEngine object
         *
         * @throw EngineException if failed to start MATLAB
         */
        std::unique_ptr<MATLABEngine> startMATLAB(const std::vector<String>& options = std::vector<String>());

        /**
        * Start MATLAB asynchronously with startup options
        *
        * @param options -  Startup options used to launch MATLAB
        * @return FutureResult<std::unique_ptr<MATLABEngine>> - A future of a unique pointer of a MATLABEngine object
        *
        * @throw none
        */
        FutureResult<std::unique_ptr<MATLABEngine>> startMATLABAsync(const std::vector<String>& options = std::vector<String>());

        /**
        * Find all shared MATLAB sessions synchronously
        *
        * @return std::vector<String> - A vector of names of shared MATLAB sessions
        * 
        * @throw none
        */
        std::vector<String> findMATLAB();

        /**
        * Find all shared MATLAB sessions synchronously
        *
        * @return FutureResult<std::vector<String>> - A future of a vector of names of shared MATLAB sessions
        *
        * @throw none
        */
        FutureResult<std::vector<String>> findMATLABAsync();

        /**
        * Start or connect to a shared MATLAB session synchronously
        *
        * @return std::unique_ptr<MATLABEngine> - A unique pointer of a MATLABEngine object
        *
        * @throw EngineException if failed to start or connect to a shared MATLAB session
        */
        std::unique_ptr<MATLABEngine> connectMATLAB();

        /**
        * Connect to a shared MATLAB session synchronously
        *
        * @param String - The name of a shared MATLAB session
        * @return std::unique_ptr<MATLABEngine> - A unique pointer of a MATLABEngine object
        *
        * @throw EngineException if failed to connect to a shared MATLAB session
        */
        std::unique_ptr<MATLABEngine> connectMATLAB(const String& name);

        /**
        * Start or connect to a shared MATLAB session asynchronously
        *
        * @return FutureResult<std::unique_ptr<MATLABEngine>> - A future of a unique pointer of a MATLABEngine object
        *
        * @throw none
        */
        FutureResult<std::unique_ptr<MATLABEngine>> connectMATLABAsync();

        /**
        * Connect to a shared MATLAB session asynchronously
        *
        * @param String - The name of a shared MATLAB session
        * @return FutureResult<std::unique_ptr<MATLABEngine>> - A future of a unique pointer of a MATLABEngine object
        *
        * @throw none
        */
        FutureResult<std::unique_ptr<MATLABEngine>> connectMATLABAsync(const String& name);

        /**
        * Terminate the engine client.  The Engine application will be no longer able to start MATLAB after calling this function.
        *
        * @return void
        *
        * @throw none
        */
        void terminateEngineClient();
    }
}


#endif  //ENGINE_FACTORY_HPP
