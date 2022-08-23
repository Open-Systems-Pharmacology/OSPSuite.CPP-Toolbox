@echo off

mkdir BuildSuiteSparse
mkdir BuildSuiteSparse\include
mkdir BuildSuiteSparse\lib64

cd src\SuiteSparse\SuiteSparse_config\
make static
rem cp -p *.a ..\..\..\BuildSuiteSparse\lib64
rem cp -p SuiteSparse_config.h ..\..\..\BuildSuiteSparse\include\

rem for SuiteSparseSubdir in AMD BTF COLAMD KLU
rem do
rem     cd ..\${SuiteSparseSubdir}\
rem     make static
rem 	cp -p Lib\*.a ..\..\..\BuildSuiteSparse\lib64
rem 	cp -p Include\*.h ..\..\..\BuildSuiteSparse\include\
rem done
rem cd ..\..\..

rem for BuildType in Debug Release
rem do
rem     cmake -BBuildSundials\${BuildType}\x64\ -Hsrc\Sundials\ -DCMAKE_BUILD_TYPE=${BuildType} -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=BuildSuiteSparse\include\ -DKLU_LIBRARY_DIR=BuildSuiteSparse\lib64\
rem     make -C BuildSundials\${BuildType}\x64\
rem     cmake -BBuild\${BuildType}\x64\ -Hsrc\OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=${BuildType} -DCMAKE_POSITION_INDEPENDENT_CODE=ON
rem     make -C Build\${BuildType}\x64\
rem     mv Build\${BuildType}\x64\libcppToolbox.so Build\${BuildType}\x64\cppToolbox.mexa64
rem done

rem nuget pack src\OSPSuite.CPP-Toolbox\OSPSuite.CPP-Toolbox_$1.nuspec -version $2
