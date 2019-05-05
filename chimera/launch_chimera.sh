#!/bin/bash
set -x

CURRENT_PATH=$1

if [[ -n $CURRENT_PATH ]]; then
    cd $CURRENT_PATH
fi
echo "Launching clang-Chimera..."
rm -Rf output
clang-chimera                                               \
        -debug -v                                           \
        -fun-op conf.csv -generate-mutants                  \
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

cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/BC12.cpp/2/axdct_report.csv > ./output/mutants/BC12.cpp/2/BC12_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/CB11.cpp/2/axdct_report.csv > ./output/mutants/CB11.cpp/2/CB11_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/BAS08.cpp/2/axdct_report.csv > ./output/mutants/BAS08.cpp/2/BAS08_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/BAS09.cpp/2/axdct_report.csv > ./output/mutants/BAS09.cpp/2/BAS09_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/BAS11.cpp/2/axdct_report.csv > ./output/mutants/BAS11.cpp/2/BAS11_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/PEA12.cpp/2/axdct_report.csv > ./output/mutants/PEA12.cpp/2/PEA12_report.csv
cat ./output/mutants/AxDCT_algorithm.cpp/2/axdct_report.csv ./output/mutants/PEA14.cpp/2/axdct_report.csv > ./output/mutants/PEA14.cpp/2/PEA14_report.csv