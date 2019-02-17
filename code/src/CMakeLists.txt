cmake_minimum_required(VERSION 3.3)
project( AxDCT )

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_subdirectory(utils)
add_subdirectory(core)
add_subdirectory(algorithms)

add_executable( displayImg main.cpp )

target_include_directories(displayImg PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(displayImg PRIVATE ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(displayImg PRIVATE ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(displayImg PRIVATE ${CMAKE_SOURCE_DIR}/include/utils)

target_link_libraries( displayImg utils ${OpenCV_LIBS} axdctCore axdctAlgorithms)