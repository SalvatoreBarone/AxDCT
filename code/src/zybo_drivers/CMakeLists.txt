find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(OpenCV_STATIC ON)

add_subdirectory(zybo_metrics)

add_library(zyboDrivers STATIC 
    bc12_zybo.cpp
    pea12_zybo.cpp
    pea14_zybo.cpp
    cb11_zybo.cpp
    bas08_zybo.cpp
    bas09_zybo.cpp
    bas11_zybo.cpp

    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/bc12_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/pea12_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/pea14_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/cb11_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/bas08_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/bas09_zybo.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/bas11_zybo.h
    
)

target_include_directories(zyboDrivers PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/zybo_drivers ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics)


add_executable( hw_eval hw_eval.cpp )

target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/mutants_eval)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/metrics)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/zybo_drivers)
target_include_directories(hw_eval PUBLIC ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics)

target_link_libraries( hw_eval
    utils 
    ${OpenCV_LIBS} 
    axdctAlgorithms
    axdctCore 
    metrics
    zyboMetrics
    zyboDrivers
)

add_executable( hw_show hw_show.cpp )

target_include_directories(hw_show PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(hw_show PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_include_directories(hw_show PUBLIC ${CMAKE_SOURCE_DIR}/include/algorithms)
target_include_directories(hw_show PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_include_directories(hw_show PUBLIC ${CMAKE_SOURCE_DIR}/include/zybo_drivers)

target_link_libraries( hw_show utils ${OpenCV_LIBS} axdctCore axdctAlgorithms zyboDrivers)
