cmake_minimum_required(VERSION 3.3)
project(AxDCT)

set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

find_package( OpenCV REQUIRED )

#if (NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()
#
#if (CMAKE_BUILD_TYPE STREQUAL "Release")
#    
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
#
#else()
#    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -ggdb3 -O0")
#
#    endif()
#    
#endif() 

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#We need the -fno-use-cxa-atexit flag for the symbol ___cxa_atexit
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -fno-use-cxa-atexit -fPIC" )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#We are going to collect all the target modules in libs folders
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
   ${PROJECT_BINARY_DIR}/libs.o/
   CACHE PATH
   "Single Directory for all dynamic and shared libraries"
   )
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/libs.a/
  CACHE PATH
  "Single Directory for all static libraries."
  )
#Then we list all of them into a the library_database file; this target must be set as dependency of the ${CMAKE_PROJECT_NAME}
add_custom_target(library_database
  COMMAND ${CMAKE_COMMAND} -E echo adding compile modules to the library_database.csv file
  COMMAND echo "\\\\\\\\list of all objects which must be loaded to run the software as JIT" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND echo "\\\\\\\\dynamic and shared libraries" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND if [ -d ${PROJECT_BINARY_DIR}/libs.o ]; then ls ${PROJECT_BINARY_DIR}/libs.o | grep 'o\\|so' >> ${PROJECT_BINARY_DIR}/library_database.csv\; fi
  COMMAND echo "\\\\\\\\static libraries" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND if [ -d ${PROJECT_BINARY_DIR}/libs.a ]; then ls ${PROJECT_BINARY_DIR}/libs.a | grep 'a' >> ${PROJECT_BINARY_DIR}/library_database.csv\; fi
)

add_subdirectory(src)

#add_dependencies(${CMAKE_PROJECT_NAME} library_database)

