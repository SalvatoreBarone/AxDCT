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

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "main.h"

using namespace cv;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    Mat originalImage = imread( argv[1], IMREAD_GRAYSCALE );
    Mat image = imread( argv[1], IMREAD_GRAYSCALE );

    assert( originalImage.data && "No image data");

    namedWindow("Original Image", WINDOW_AUTOSIZE );
    imshow("Original Image", originalImage);

    //TODO: In image ho la matrice dei pixel dell'immagine

    namedWindow("Modified Image", WINDOW_AUTOSIZE );
    imshow("Modified Image", image);
    
    waitKey(0);
    return 0;
}