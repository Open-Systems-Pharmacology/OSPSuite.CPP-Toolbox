cmake_minimum_required (VERSION 3.9)
project (OSPSuite.CPPToolbox)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(LibXml2 REQUIRED)

include_directories (
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/include 
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../Sundials/include
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/extern/include
	${LIBXML2_INCLUDE_DIR}
)

file (GLOB SOURCES 
        ${OSPSuite.CPPToolbox_SOURCE_DIR}/src/*.cpp
)

add_library (OSPSuite.CPPToolbox SHARED ${SOURCES})

target_link_libraries (OSPSuite.CPPToolbox 
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmat.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmex.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/bin/glnxa64/libmx.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/runtime/glnxa64/libmwmclmcrrt.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../lib/Matlab2017b_Linux/runtime/glnxa64/libmwmclmcrrt.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_arkode.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_cvode.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_cvodes.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_ida.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_idas.so
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../Build/${CMAKE_BUILD_TYPE}/x64/libsundials_kinsol.so
	${LIBXML2_LIBRARIES}
)