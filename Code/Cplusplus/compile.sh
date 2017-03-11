#!/bin/bash

# enable c++ 11 on centos 6
# scl enable devtoolset-2 bash

# Set compilation mode to use
RELEASE=0
DEBUG=1
COMPILATION=$RELEASE
#COMPILATION=$DEBUG

GIT_DIRECTORY=`git rev-parse --show-toplevel`

# Set compilation flags
if [ "$COMPILATION" -eq "$RELEASE" ] ; then
    FLAGS="-O3 -DNDEBUG -m64 -std=c++11"
elif [ "$COMPILATION" -eq "$DEBUG" ] ; then
    FLAGS="-O0 -g -m64 -std=c++11"
else
    echo "Unknown compilation mode"
fi

# Using g++
g++ -I${GIT_DIRECTORY}/Code/Cplusplus -I${GIT_DIRECTORY}/Code/Cplusplus/Eigen327 ${FLAGS} -fPIC -c -o HeteroGauss.o heterogauss.cpp

# Make a library
ar rcs libHeteroGauss.a HeteroGauss.o


