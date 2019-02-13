#!/bin/bash
set -x

echo "Launching clang-Chimera..."
clang-chimera -debug -v -fun-op conf.csv -generate-mutants ../code/src/main.cpp -o output -- -std=c++11 -I../code/include -I/usr/local/include/opencv4 -lopencv_core
