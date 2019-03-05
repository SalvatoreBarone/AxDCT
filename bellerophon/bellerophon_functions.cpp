/******************************************************************************
 * @file   bellerophon_functions.cpp
 * @author Andrea Aletto
 * @date   13 Feb 2019
 * @brief  Bellerophon functions for the DCT Algorithm Approximated through
 *         clang-chimera
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>


#include "algorithms_list.h"
#include "algorithms/BAS08.h"
#include "core/dct.h"
#include "core/metrics.h"


char oracle_path[]= "./main_oracle.txt";

// extern "C"
// {
//     void *__dso_handle = NULL;
// }

extern "C" double BELLERO_getError()
{

    const std::string img_path = "/home/andrea/lena512color.bmp";
    const cv::Mat bgrImg = imread( img_path.c_str(), cv::IMREAD_COLOR );

    cv::Mat dst = bgrImg;
    // Direct and inverse transform
    transformImage(bgrImg,dst, new BC12 );
    inverseTransformImage(dst, dst, new BC12);
    

    // Compute PSNR and return it
    return compute_mse(bgrImg, dst);

}


