cmake_minimum_required(VERSION 3.3)

add_library(axdctCore STATIC 
    
    dct.cpp
    metrics.cpp
    ${CMAKE_SOURCE_DIR}/include/utils/mat_operations.h
    ${CMAKE_SOURCE_DIR}/include/core/dct.h
    ${CMAKE_SOURCE_DIR}/include/core/metrics.h
)

target_include_directories(axdctCore PUBLIC ${CMAKE_SOURCE_DIR}/include/core ${CMAKE_SOURCE_DIR}/include/utils)
#set_target_properties(axdctCore PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/include/core/dct.h")

target_link_libraries(axdctCore PUBLIC ${OpenCV_LIBS} utils)