//
// AxDCT - A collection of inexact DCT Algorithms
// Copyright (C) 2019 Andrea Aletto <andrea.aletto8@gmail.com>
//
// This file is part of AxDCT.
//
// AxDCT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// AxDCT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with AxDCT.  If not, see <http://www.gnu.org/licenses/>.
//

/******************************************************************************
 * @file   main.cpp
 * @author Andrea Aletto
 * @date   4 feb 2019
 * @brief  Implementation of main executable functions
 ******************************************************************************/

#include "core/dct.h"
#include "algorithms_list.h"

#define CHECKPOINT (std::cerr<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    // Load img
    cv::Mat bgrImg = imread( argv[1], cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");

    // Declare an empty image for transformation
    cv::Mat transfImg = bgrImg;

    // Direct and inverse transform
    transformImage(bgrImg,transfImg, new BC12 );
    inverseTransformImage(transfImg, bgrImg, new BC12);

    // Show the approximate image 
    cv::namedWindow("Approximate Image", cv::WINDOW_AUTOSIZE );
    imshow("Approximate Image", bgrImg);


    // Load img
    cv::Mat bgrImg2 = imread( argv[1], cv::IMREAD_COLOR );
    assert( bgrImg2.data && "No image data");

    // Declare an empty image for transformation
    cv::Mat transfImg2 = bgrImg2;

    // Direct and inverse transform
    transformImage(bgrImg2,transfImg2);
    inverseTransformImage(transfImg2, bgrImg2);

    // Show the approximate image 
    cv::namedWindow("Approximate Image 2", cv::WINDOW_AUTOSIZE );
    imshow("Approximate Image 2", bgrImg2);
    
    cv::waitKey(0);
    return 0;
}








