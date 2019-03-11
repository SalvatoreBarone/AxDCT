find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(OpenCV_STATIC ON)

add_library(axdctMutatedAlgorithms OBJECT 
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/AxDCT_algorithm.cpp/1/AxDCT_algorithm.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BC12.cpp/1/BC12.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/CB11.cpp/1/CB11.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS09.cpp/1/BAS09.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS08.cpp/1/BAS08.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS11.cpp/1/BAS11.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/PEA12.cpp/1/PEA12.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/PEA14.cpp/1/PEA14.cpp

    ${CMAKE_SOURCE_DIR}/include/mutants_eval/nablist.h
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

target_include_directories(axdctMutatedAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms ${CMAKE_SOURCE_DIR}/include)

add_executable( psnr_eval psnr_eval.cpp )

target_include_directories(psnr_eval PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(psnr_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(psnr_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(psnr_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_include_directories(psnr_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/mutants_eval)

target_link_libraries( psnr_eval
    utils 
    ${OpenCV_LIBS} 
    axdctCore 
    axdctMutatedAlgorithms
    /usr/lib/libInexactAdders.a
)

add_dependencies(psnr_eval library_database)


#Then we list all of them into a the library_database file
add_custom_target(library_database ALL
  COMMAND ${CMAKE_COMMAND} -E echo adding compile modules to the library_database.csv file
  COMMAND echo "\\\\\\\\list of all objects which must be loaded to run the software as JIT" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND echo "\\\\\\\\dynamic and shared libraries" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND if [ -d ${PROJECT_BINARY_DIR}/libs.o ]; then ls ${PROJECT_BINARY_DIR}/libs.o | grep 'o\\|so' >> ${PROJECT_BINARY_DIR}/library_database.csv\; fi
  COMMAND echo "\\\\\\\\static libraries" >> ${PROJECT_BINARY_DIR}/library_database.csv
  COMMAND if [ -d ${PROJECT_BINARY_DIR}/libs.a ]; then ls ${PROJECT_BINARY_DIR}/libs.a | grep 'a' >> ${PROJECT_BINARY_DIR}/library_database.csv\; fi
  DEPENDS utils axdctCore axdctAlgorithms
)
