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
 * @file   dct.h
 * @author Andrea Aletto
 * @date   16 feb 2019
 * @brief  Declaration of dct functions in approximate context
 ******************************************************************************/

#ifndef _DCT_H
#define _DCT_H

#include <iostream>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include "../utils/mat_operations.h"
#include "../algorithms_list.h"
#include "../user_defines.h"

void AxDCT(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);

void  y_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);
void cb_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);
void cr_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);

void  y_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);
void cb_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);
void cr_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);

void transformImage(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);
void inverseTransformImage(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg = new __USER_DEFAULT_ALGORITHM);

#endif /* _DCT_H */
