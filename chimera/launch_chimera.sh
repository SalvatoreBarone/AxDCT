#!/bin/bash
set -x

CURRENT_PATH=$1

if [[ -n $CURRENT_PATH ]]; then
    cd $CURRENT_PATH
fi
echo "Launching clang-Chimera..."
rm -Rf output
clang-chimera -debug -v -fun-op conf.csv -generate-mutants  \
        ../code/src/algorithms/AxDCT_algorithm.cpp          \
        ../code/src/algorithms/BAS08.cpp                    \
        ../code/src/algorithms/BAS09.cpp                    \
        ../code/src/algorithms/BAS11.cpp                    \
        ../code/src/algorithms/BC12.cpp                     \
        ../code/src/algorithms/CB11.cpp                     \
        ../code/src/algorithms/PEA12.cpp                    \
        ../code/src/algorithms/PEA14.cpp                    \
        -o output --                                        \
        -std=c++11                                          \
        -I../code/include                                   \
        -I../code/include/algorithms                        \
        -I../code/include/utils                             \
        -I../code/include/core                              \
        -I/usr/local/include/opencv4                        \
        -lopencv_core   
