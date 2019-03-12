cmake_minimum_required(VERSION 3.3)

add_library(metrics STATIC 
    
    metrics.cpp
    psnr_metric_eval.cpp
    mse_metric_eval.cpp
    ad_metric_eval.cpp
    md_metric_eval.cpp
    ${CMAKE_SOURCE_DIR}/include/utils/mat_operations.h
    ${CMAKE_SOURCE_DIR}/include/core/dct.h
    ${CMAKE_SOURCE_DIR}/include/metrics/metrics.h
)

target_include_directories(metrics PUBLIC ${CMAKE_SOURCE_DIR}/include/core ${CMAKE_SOURCE_DIR}/include/utils ${CMAKE_SOURCE_DIR}/include/metrics)

target_link_libraries(metrics PUBLIC ${OpenCV_LIBS} utils axdctCore)