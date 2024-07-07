#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Missing or invalid arguments supplied"
    echo "Usage: createNugetPackageLinux.sh <version>"
    echo " e.g.: createNugetPackageLinux.sh 1.0.1"
    exit 1
fi

nuget pack OSPSuite.CPP-Toolbox_Linux.nuspec -version $1
