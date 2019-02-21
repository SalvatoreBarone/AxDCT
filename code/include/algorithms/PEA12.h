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
 * @file   PEA12.h
 * @author Andrea Aletto
 * @date   21 feb 2019
 * @brief  Declaration of PEA12 algorithm class
 ******************************************************************************/

#ifndef _PEA12_H
#define _PEA12_H

#include <opencv2/opencv.hpp>
#include "AxDCT_algorithm.h"

class PEA12 : public AxDCT_algorithm
{

    private:
        virtual void dct1d(const cv::Mat& input, cv::Mat& output);

        virtual cv::Mat getYQuantizationMatix();
        virtual cv::Mat getCrQuantizationMatix();
        virtual cv::Mat getCbQuantizationMatix();

        virtual cv::Mat getYDequantizationMatix();
        virtual cv::Mat getCrDequantizationMatix();
        virtual cv::Mat getCbDequantizationMatix();

        cv::Mat getT();
        cv::Mat getD();
        cv::Mat getQ();
        cv::Mat getCQ();

    public:
        PEA12() : AxDCT_algorithm() {};

};

#endif /* _PEA12_H */