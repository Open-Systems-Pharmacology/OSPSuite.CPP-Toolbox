#! /bin/sh

#call: buildLinux.sh distributionName version
# e.g. buildLinux.sh CentOS7 4.0.0.49

cd src/SuiteSparse/
#make clean
make library
cd ../..

cmake -BBuildSundials/Release/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Release -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=src/SuiteSparse/include -DKLU_LIBRARY_DIR=src/SuiteSparse/lib
make -C BuildSundials/Release/x64/
cmake -BBuild/Release/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
make -C Build/Release/x64/
cp -p Build/Release/x64/libOSPSuite.CPP-Toolbox.so Build/Release/x64/libOSPSuite.CPP-Toolbox.mexa64

cmake -BBuildSundials/Debug/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Debug -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DKLU_INCLUDE_DIR=src/SuiteSparse/include -DKLU_LIBRARY_DIR=src/SuiteSparse/lib
make -C BuildSundials/Debug/x64/
cmake -BBuild/Debug/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=Debug -DCMAKE_POSITION_INDEPENDENT_CODE=ON
make -C Build/Debug/x64/
cp -p Build/Debug/x64/libOSPSuite.CPP-Toolbox.so Build/Debug/x64/libOSPSuite.CPP-Toolbox.mexa64

nuget pack src/OSPSuite.CPP-Toolbox/OSPSuite.CPP-Toolbox_$1.nuspec -version $2
