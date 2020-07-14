cmake_minimum_required (VERSION 3.9)
project (OSPSuite.CPPToolbox)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_SHARED_LINKER_FLAGS "-static-libstdc++ -Wl,-Bdynamic") #cppToolbox 15
#set (CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++ -Wl,-Bdynamic") #cppToolbox 12
#set (CMAKE_SHARED_LINKER_FLAGS "-Wl,-Bdynamic") #cppToolbox 11
#set (CMAKE_SHARED_LINKER_FLAGS "-lrt -lm -ldl -static-libgcc -static-libstdc++ -Wl,-Bdynamic") #cppToolbox 10
#set (CMAKE_SHARED_LINKER_FLAGS "-lrt -lm -ldl -static-libgcc -static-libstdc++ -Wl,-Bdynamic")
#set (CMAKE_LINKER_FLAGS "-lrt -lm -ldl")
#set (CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

include_directories (
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/include 
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../Sundials/include
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/extern/include
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/include
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/include
)

file (GLOB SOURCES 
        ${OSPSuite.CPPToolbox_SOURCE_DIR}/src/*.cpp
)

add_library (cppToolbox SHARED ${SOURCES})

target_link_libraries (cppToolbox 
	dl rt m pthread
	c gcc_s gcc 
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmat.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmex.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmx.so
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/runtime/glnxa64/libmwmclmcrrt.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/arkode/libsundials_arkode.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/cvode/libsundials_cvode.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/cvodes/libsundials_cvodes.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/ida/libsundials_ida.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/idas/libsundials_idas.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/kinsol/libsundials_kinsol.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/nvector/manyvector/libsundials_nvecmanyvector.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/nvector/serial/libsundials_nvecserial.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/band/libsundials_sunlinsolband.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/dense/libsundials_sunlinsoldense.a	
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/klu/libsundials_sunlinsolklu.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/pcg/libsundials_sunlinsolpcg.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/spbcgs/libsundials_sunlinsolspbcgs.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/spfgmr/libsundials_sunlinsolspfgmr.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/spgmr/libsundials_sunlinsolspgmr.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunmatrix/band/libsundials_sunmatrixband.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunmatrix/dense/libsundials_sunmatrixdense.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunnonlinsol/fixedpoint/libsundials_sunnonlinsolfixedpoint.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/AMD/Lib/libamd.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/BTF/Lib/libbtf.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/COLAMD/Lib/libcolamd.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/KLU/Lib/libklu.a
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
)