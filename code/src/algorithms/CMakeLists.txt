cmake_minimum_required(VERSION 3.3)


add_library(axdctAlgorithms OBJECT 
    AxDCT_algorithm.cpp
    BC12.cpp
    CB11.cpp
    BAS09.cpp
    BAS08.cpp
    BAS11.cpp
    PEA12.cpp
    PEA14.cpp

    ${CMAKE_SOURCE_DIR}/include/user_defines.h
    ${CMAKE_SOURCE_DIR}/include/algorithms_list.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/AxDCT_algorithm.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/BC12.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/CB11.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/BAS09.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/BAS08.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/BAS11.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/PEA12.h
    ${CMAKE_SOURCE_DIR}/include/algorithms/PEA14.h
    
)

target_include_directories(axdctAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms ${CMAKE_SOURCE_DIR}/include)