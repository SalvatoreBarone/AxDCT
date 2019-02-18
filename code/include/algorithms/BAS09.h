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
 * @file   BAS09.h
 * @author Andrea Aletto
 * @date   11 feb 2019
 * @brief  Declaration of BAS09 algorithm class
 ******************************************************************************/

#ifndef _BAS09_H
#define _BAS09_H

#include <opencv2/opencv.hpp>
#include "AxDCT_algorithm.h"

class BAS09 : public AxDCT_algorithm
{

    private:
        void dct1d(const cv::Mat& input, cv::Mat& output);

        cv::Mat getT();
        cv::Mat getD();
        cv::Mat getQ();
        cv::Mat getCQ();

    public:
        BAS09() : AxDCT_algorithm() {};
        void y_quantizate(const cv::Mat&, cv::Mat&);
        void cr_quantizate(const cv::Mat&, cv::Mat&);
        void cb_quantizate(const cv::Mat&, cv::Mat&);

        void y_dequantizate(const cv::Mat&, cv::Mat&);
        void cr_dequantizate(const cv::Mat&, cv::Mat&);
        void cb_dequantizate(const cv::Mat&, cv::Mat&);
};

#endif /* _BAS09_H */