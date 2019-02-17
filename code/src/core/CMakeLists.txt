cmake_minimum_required(VERSION 3.3)

add_library(axdctCore SHARED 
    
    dct.cpp
    ${CMAKE_SOURCE_DIR}/include/utils/mat_operations.h
    ${CMAKE_SOURCE_DIR}/include/core/dct.h
)

target_include_directories(axdctCore PRIVATE ${CMAKE_SOURCE_DIR}/include/core ${CMAKE_SOURCE_DIR}/include/utils)
set_target_properties(axdctCore PROPERTIES PRIVATE_HEADER "${CMAKE_SOURCE_DIR}/include/core/dct.h")

target_link_libraries(axdctCore PRIVATE ${OpenCV_LIBS} utils)