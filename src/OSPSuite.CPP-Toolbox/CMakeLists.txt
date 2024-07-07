cmake_minimum_required (VERSION 3.9)
project (OSPSuite.CPPToolbox)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
if (UNIX)
set (CMAKE_SHARED_LINKER_FLAGS "-static-libstdc++ -Wl,-Bdynamic") #-static-libgcc
endif (UNIX)

if (UNIX)
include_directories (
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/include 
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Ubuntu22/CVODES/include
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSuiteSparse/include
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/include
)

link_directories (
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSuiteSparse/lib64
${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Ubuntu22/CVODES/bin/native/x64/${CMAKE_BUILD_TYPE}
)

endif (UNIX)

if (WIN32)
include_directories (
    ${OSPSuite.CPPToolbox_SOURCE_DIR}/include 
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Windows/CVODES/include
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSuiteSparse/include
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSundials/${CMAKE_BUILD_TYPE}/x64/include
)

link_directories (
#	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../BuildSuiteSparse/lib64
${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Windows/CVODES/bin/native/x64/${CMAKE_BUILD_TYPE}
)

endif (WIN32)


file (GLOB SOURCES 
        ${OSPSuite.CPPToolbox_SOURCE_DIR}/src/*.cpp
)

add_library (cppToolbox SHARED ${SOURCES})

if (UNIX)
file(GLOB CVODE_LIBS "${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Ubuntu22/CVODES/bin/native/x64/${CMAKE_BUILD_TYPE}/*.a")

target_link_libraries (cppToolbox 
    dl rt m pthread
	c gcc_s gcc 
#	klu amd btf colamd suitesparseconfig
#	suitesparseconfig amd btf colamd klu
${CVODE_LIBS}
# ${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Ubuntu22/CVODES/bin/native/x64/${CMAKE_BUILD_TYPE}/*.a
)
endif (UNIX)

if (WIN32)
#target_link_options(cppToolbox PRIVATE
#        -sEXPORTED_FUNCTIONS=['mexFunction'])
        
target_link_libraries (cppToolbox 
#	klu amd btf colamd suitesparseconfig
	${OSPSuite.CPPToolbox_SOURCE_DIR}/../../packages/CVODES.Windows/CVODES/bin/native/x64/${CMAKE_BUILD_TYPE}/*.lib
)
endif (WIN32)