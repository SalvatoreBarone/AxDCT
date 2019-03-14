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
 * @file   metrics.h
 * @author Andrea Aletto
 * @date   28 feb 2019
 * @brief  Declaration of evaluating metrics for transformed images
 ******************************************************************************/

#ifndef _METRICS_H_
#define _METRICS_H_

#include <iostream>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include "../utils/mat_operations.h"
#include "../algorithms_list.h"
#include "../user_defines.h"

// double compute_mse(const cv::Mat& orig, const cv::Mat& target, int component);
// double compute_psnr(const cv::Mat& orig, const cv::Mat& target, int component);


double compute_mse(const cv::Mat& orig, const cv::Mat& target);
double compute_reduction(const double exact_param, const double inexact_param, const int nab, const int n_bit = 16);
double compute_ad(const cv::Mat& orig, const cv::Mat& target);
double compute_md(const cv::Mat& orig, const cv::Mat& target);

double compute_psnr(const cv::Mat& orig, const cv::Mat& target);
double compute_mssim(const cv::Mat& orig, const cv::Mat& target);

#endif
