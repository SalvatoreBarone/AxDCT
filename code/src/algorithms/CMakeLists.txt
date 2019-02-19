cmake_minimum_required(VERSION 3.3)


add_library(axdctAlgorithms STATIC 
    AxDCT_algorithm.cpp
    BC12.cpp
    BAS09.cpp
    BAS08.cpp
)

target_include_directories(axdctAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)