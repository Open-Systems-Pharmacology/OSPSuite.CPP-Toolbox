#! /bin/sh

if [ "$#" -ne 2 ]; then
    echo "Missing or invalid arguments supplied"
    echo "Usage: buildLinux.sh <distributionName> <version>"
    echo " e.g.: buildLinux.sh Ubuntu18 1.0.1"
    exit 1
fi

mkdir BuildSuiteSparse
mkdir BuildSuiteSparse/include
mkdir BuildSuiteSparse/lib64

cd src/SuiteSparse/SuiteSparse_config/
make static
cp -p *.a ../../../BuildSuiteSparse/lib64
cp -p SuiteSparse_config.h ../../../BuildSuiteSparse/include/

for SuiteSparseSubdir in AMD BTF COLAMD KLU
do
    cd ../${SuiteSparseSubdir}/
    make static
	cp -p Lib/*.a ../../../BuildSuiteSparse/lib64
	cp -p Include/*.h ../../../BuildSuiteSparse/include/
done
cd ../../..

for BuildType in Debug Release
do
    cmake -BBuildSundials/${BuildType}/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=${BuildType} -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=BuildSuiteSparse/include/ -DKLU_LIBRARY_DIR=BuildSuiteSparse/lib64/
    make -C BuildSundials/${BuildType}/x64/
    cmake -BBuild/${BuildType}/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=${BuildType} -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    make -C Build/${BuildType}/x64/
    mv Build/${BuildType}/x64/libcppToolbox.so Build/${BuildType}/x64/cppToolbox.mexa64
done

nuget pack src/OSPSuite.CPP-Toolbox/OSPSuite.CPP-Toolbox_$1.nuspec -version $2
