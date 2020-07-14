#! /bin/sh

#call: buildLinux.sh <distributionName> <version>
# e.g. buildLinux.sh CentOS7 1.0.1

cd src/SuiteSparse/SuiteSparse_config/
make library
cd ../AMD/
make library
cd ../BTF/
make library
cd ../COLAMD/
make library
cd ../KLU/
make library
cd ../../..

mkdir BuildSuiteSparse
mkdir BuildSuiteSparse/include
mkdir BuildSuiteSparse/lib64

cp -p src/SuiteSparse/AMD/Lib/*.a BuildSuiteSparse/lib64
cp -p src/SuiteSparse/BTF/Lib/*.a BuildSuiteSparse/lib64
cp -p src/SuiteSparse/COLAMD/Lib/*.a BuildSuiteSparse/lib64
cp -p src/SuiteSparse/KLU/Lib/*.a BuildSuiteSparse/lib64
cp -p src/SuiteSparse/SuiteSparse_config/*.a BuildSuiteSparse/lib64

cp -p src/SuiteSparse/AMD/Include/* BuildSuiteSparse/include/
cp -p src/SuiteSparse/BTF/Include/* BuildSuiteSparse/include/
cp -p src/SuiteSparse/COLAMD/Include/* BuildSuiteSparse/include/
cp -p src/SuiteSparse/KLU/Include/* BuildSuiteSparse/include/
cp -p src/SuiteSparse/SuiteSparse_config/SuiteSparse_config.h BuildSuiteSparse/include/

for BuildType in Debug Release
do
    cmake -BBuildSundials/${BuildType}/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=${BuildType} -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=BuildSuiteSparse/include/ -DKLU_LIBRARY_DIR=BuildSuiteSparse/lib64/
    make -C BuildSundials/${BuildType}/x64/
    cmake -BBuild/${BuildType}/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=${BuildType} -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    make -C Build/${BuildType}/x64/
    cp -p Build/${BuildType}/x64/libOSPSuite.CPP-Toolbox.so Build/${BuildType}/x64/cppToolbox.mexa64
done

#cmake -BBuildSundials/Release/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Release -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=BuildSuiteSparse/include/ -DKLU_LIBRARY_DIR=BuildSuiteSparse/lib64/

#cmake -BBuildSundials/Release/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Release -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=src/SuiteSparse/include -DKLU_LIBRARY_DIR=src/SuiteSparse/lib
#make -C BuildSundials/Release/x64/
#cmake -BBuild/Release/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
#make -C Build/Release/x64/
#cp -p Build/Release/x64/libOSPSuite.CPP-Toolbox.so Build/Release/x64/libOSPSuite.CPP-Toolbox.mexa64

#cmake -BBuildSundials/Debug/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Debug -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=src/SuiteSparse/include -DKLU_LIBRARY_DIR=src/SuiteSparse/lib
#make -C BuildSundials/Debug/x64/
#cmake -BBuild/Debug/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=Debug -DCMAKE_POSITION_INDEPENDENT_CODE=ON
#make -C Build/Debug/x64/
#cp -p Build/Debug/x64/libOSPSuite.CPP-Toolbox.so Build/Debug/x64/libOSPSuite.CPP-Toolbox.mexa64

nuget pack src/OSPSuite.CPP-Toolbox/OSPSuite.CPP-Toolbox_$1.nuspec -version $2
