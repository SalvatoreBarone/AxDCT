find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(OpenCV_STATIC ON)

add_library(zyboMetrics STATIC 
    
    ${CMAKE_SOURCE_DIR}/src/metrics/metrics.cpp
    zybo_psnr_metric_eval.cpp
    zybo_mse_metric_eval.cpp
    zybo_ad_metric_eval.cpp
    zybo_md_metric_eval.cpp
    zybo_mssim_metric_eval.cpp
    zybo_dssim_metric_eval.cpp
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_psnr_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_mse_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_ad_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_md_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_mssim_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics/zybo_dssim_metric_eval.h
    ${CMAKE_SOURCE_DIR}/include/utils/mat_operations.h
    ${CMAKE_SOURCE_DIR}/include/core/dct.h
    ${CMAKE_SOURCE_DIR}/include/metrics/metrics.h
)

target_include_directories(zyboMetrics PUBLIC ${CMAKE_SOURCE_DIR}/include/core ${CMAKE_SOURCE_DIR}/include/utils ${CMAKE_SOURCE_DIR}/include/metrics ${CMAKE_SOURCE_DIR}/include/zybo_drivers ${CMAKE_SOURCE_DIR}/include/zybo_drivers/zybo_metrics ${CMAKE_SOURCE_DIR}/include/algorithms)

target_link_libraries(zyboMetrics PUBLIC ${OpenCV_LIBS} utils axdctCore axdctAlgorithms)