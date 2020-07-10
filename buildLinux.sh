#! /bin/sh

#call: buildLinux.sh distributionName version
# e.g. buildLinux.sh CentOS7 4.0.0.49

#git submodule update --init --recursive

cd src/SuiteSparse/
#make clean
make static
cd ../..
cmake -BBuildSundials/Release/x64/ -Hsrc/Sundials/ -DCMAKE_BUILD_TYPE=Release -DEXAMPLES_ENABLE_C=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DKLU_ENABLE=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON
make -C BuildSundials/Release/x64/
cmake -BBuild/Release/x64/ -Hsrc/CPP-Toolbox -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON
make -C Build/Release/x64/


# cp -p Build/Release/x64/src/*/*.a Build/Release/x64/

#nuget pack src/OSPSuite.SimModelSolver_CVODES/OSPSuite.SimModelSolver_CVODES_$1.nuspec -version $2
