#!/bin/bash
set -x

echo "Launching clang-Chimera..."
clang-chimera -debug -v -fun-op conf.csv -generate-mutants ../code/src/algorithms/BC12.cpp -o output -- -std=c++11 -I../code/include -I../code/include/algorithms -I../code/include/utils -I../code/include/core -I/usr/local/include/opencv4 -lopencv_core
