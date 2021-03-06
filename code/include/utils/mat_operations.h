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
 * @file   mat_operations.h
 * @author Andrea Aletto
 * @date   16 feb 2019
 * @brief  Declaration of utils functions for image structures handling
 ******************************************************************************/

#ifndef _MAT_OPERATIONS_H
#define _MAT_OPERATIONS_H

#include <iostream>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>

template<typename T>
void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type = CV_64FC1);
// template<typename T>
// void matrix_mult(cv::Mat const&, cv::Mat const&, cv::Mat&, int);

cv::Mat **splitInTiles(const cv::Mat &input, int blockSize);

cv::Mat mergeTiles( cv::Mat **tiles, int imgWidth, int imgLength, int blockSize = 8, bool deallocTiles = true);


extern template void matrix_mult<unsigned char>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);
extern template void matrix_mult<int16_t>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);
extern template void matrix_mult<double>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);

#endif /* _MAT_OPERATIONS_H */
