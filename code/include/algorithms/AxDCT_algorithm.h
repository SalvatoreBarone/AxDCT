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
 * @file   AxDCT_algorithm.h
 * @author Andrea Aletto
 * @date   18 feb 2019
 * @brief  AxDCT algorithm base class declaration
 ******************************************************************************/

#ifndef _AXDCT_ALGORITHM_H
#define _AXDCT_ALGORITHM_H

#include <iostream>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <inexact_adders.h>
#include "../utils/mat_operations.h"

class AxDCT_algorithm{
    public:
        virtual void dct(const cv::Mat&, cv::Mat&);
        virtual void y_quantizate(const cv::Mat&, cv::Mat&) = 0;
        virtual void cr_quantizate(const cv::Mat&, cv::Mat&) = 0;
        virtual void cb_quantizate(const cv::Mat&, cv::Mat&) = 0;

        virtual void y_dequantizate(const cv::Mat&, cv::Mat&) = 0;
        virtual void cr_dequantizate(const cv::Mat&, cv::Mat&) = 0;
        virtual void cb_dequantizate(const cv::Mat&, cv::Mat&) = 0;

    protected:
        virtual void dct1d(const cv::Mat&, cv::Mat&) = 0;

};

#endif /* _AXDCT_ALGORITHM_H */