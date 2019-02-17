cmake_minimum_required(VERSION 3.3)


add_library(axdctAlgorithms SHARED 
    BC12.cpp
    BAS09.cpp
)

target_include_directories(axdctAlgorithms PRIVATE ${CMAKE_SOURCE_DIR}/include/algorithms)