cmake_minimum_required(VERSION 3.3)
project(AxDCT)

set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package( OpenCV REQUIRED )

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fno-use-cxa-atexit -fPIC")

else()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -ggdb3 -O0 -fno-use-cxa-atexit -fPIC")

    endif()
    
endif() 

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs.o )
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs.a )

add_subdirectory(src)


