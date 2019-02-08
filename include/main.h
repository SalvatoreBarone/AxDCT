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
 * @file   main.h
 * @author Andrea Aletto
 * @date   4 feb 2019
 * @brief  Declaration of main executable functions
 ******************************************************************************/
#ifndef _MAIN_H
#define _MAIN_H

#include <iostream>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>

enum AxDCT_Algorithm{
    BAS08,
    BAS09,
    BC12
};

void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type = CV_8U);
cv::Mat **splitInTiles(const cv::Mat &input, int blockSize);
cv::Mat mergeTiles( cv::Mat **tiles, int imgWidth, int imgLength, int blockSize = 8, bool deallocTiles = true);
void AxDCT(const cv::Mat& tile, const cv::Mat& T, cv::Mat& output);
void retrieveParameters(const AxDCT_Algorithm alg, cv::Mat& T);



#endif