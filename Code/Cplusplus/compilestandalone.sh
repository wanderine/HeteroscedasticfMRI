#!/bin/bash

GIT_DIRECTORY=`git rev-parse --show-toplevel`

# Set compilation mode to use
RELEASE=0
DEBUG=1
COMPILATION=$RELEASE
#COMPILATION=$DEBUG

# Set compilation flags
if [ "$COMPILATION" -eq "$RELEASE" ] ; then
    FLAGS="-O3 -DNDEBUG -m64 -std=c++11 -fopenmp"
elif [ "$COMPILATION" -eq "$DEBUG" ] ; then
    FLAGS="-O0 -g -m64 -std=c++11"
else
    echo "Unknown compilation mode"
fi

g++ HeteroGLM.cpp -I${GIT_DIRECTORY}/Code/Cplusplus -I${GIT_DIRECTORY}/Code/Cplusplus/Eigen327  -L${GIT_DIRECTORY}/Code/Cplusplus -L${GIT_DIRECTORY}/Code/Cplusplus/nifticlib-2.0.0/lib -I${GIT_DIRECTORY}/Code/Cplusplus/nifticlib-2.0.0/niftilib -I${GIT_DIRECTORY}/Code/Cplusplus/nifticlib-2.0.0/znzlib ${FLAGS}  -lHeteroGauss -lniftiio -lznz -lz  -o HeteroGLM
