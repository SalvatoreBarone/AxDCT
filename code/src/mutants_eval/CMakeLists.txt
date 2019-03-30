find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(OpenCV_STATIC ON)

add_library(axdctMutatedAlgorithms OBJECT 
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/AxDCT_algorithm.cpp/2/AxDCT_algorithm.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BC12.cpp/2/BC12.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/CB11.cpp/2/CB11.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS09.cpp/2/BAS09.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS08.cpp/2/BAS08.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/BAS11.cpp/2/BAS11.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/PEA12.cpp/2/PEA12.cpp
    ${CMAKE_SOURCE_DIR}/../chimera/output/mutants/PEA14.cpp/2/PEA14.cpp

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

add_executable( mutant_eval mutant_eval.cpp )

target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/mutants_eval)
target_include_directories(mutant_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/metrics)

target_link_libraries( mutant_eval
    utils 
    ${OpenCV_LIBS} 
    axdctCore 
    axdctMutatedAlgorithms
    metrics
    /usr/lib/libInexactAdders.a
)

add_executable( mutant_show mutant_show.cpp )

target_include_directories(mutant_show PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(mutant_show PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(mutant_show PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(mutant_show PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_include_directories(mutant_show PUBLIC ${CMAKE_SOURCE_DIR}/include/mutants_eval)

target_link_libraries( mutant_show
    utils 
    ${OpenCV_LIBS} 
    axdctCore 
    axdctMutatedAlgorithms
    /usr/lib/libInexactAdders.a
)

