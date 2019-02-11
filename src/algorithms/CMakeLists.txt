cmake_minimum_required(VERSION 3.3)


add_library(axdctAlgorithms STATIC 
    BC12.cpp
)

target_include_directories(axdctAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)