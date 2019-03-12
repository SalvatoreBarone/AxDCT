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
 * @file   mse_metric_eval.h
 * @author Andrea Aletto
 * @date   12 mar 2019
 * @brief  Declaration of mse calc function
 ******************************************************************************/

#ifndef _MSE_METRIC_EVAL_H
#define _MSE_METRIC_EVAL_H

#include "../core/dct.h"
#include "../algorithms_list.h"
#include "metrics.h"

namespace metrics {
    double BC12_MSE(const cv::Mat& orig);
    double CB11_MSE(const cv::Mat& orig);
    double BAS08_MSE(const cv::Mat& orig);
    double BAS09_MSE(const cv::Mat& orig);
    double BAS11_MSE(const cv::Mat& orig, double a_param);
    double PEA12_MSE(const cv::Mat& orig);
    double PEA14_MSE(const cv::Mat& orig);
}

#endif /* _MSE_METRIC_EVAL_H */