#! /bin/sh

for BuildType in Debug Release
do
    cmake -BBuildLinux/${BuildType}/x64/ -Hsrc/OSPSuite.CPP-Toolbox -DCMAKE_BUILD_TYPE=${BuildType} -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    make -C BuildLinux/${BuildType}/x64/
    
    mkdir -p Dist/Linux/${BuildType}/x64
    cp -p BuildLinux/${BuildType}/x64/*.so Dist/Linux/${BuildType}/x64/
done

