cmake_minimum_required(VERSION 3.3)

add_library(utils STATIC 
    
    mat_operations.cpp
    ${CMAKE_SOURCE_DIR}/include/utils/mat_operations.h
)

target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)

target_link_libraries(utils PUBLIC ${OpenCV_LIBS} )