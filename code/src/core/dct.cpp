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
 * @file   dct.cpp
 * @author Andrea Aletto
 * @date   16 feb 2019
 * @brief  Implementation of dct functions in approximate context
 ******************************************************************************/

#include "dct.h"

void AxDCT(const cv::Mat& tile, cv::Mat& output){

    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.dct(tile, output);
    delete &alg;

}

void y_quantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.y_quantizate(tile, output);
    delete &alg;
}

void cr_quantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.cr_quantizate(tile, output);
    delete &alg;
}

void cb_quantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.cb_quantizate(tile, output);
    delete &alg;
}

void y_dequantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.y_dequantizate(tile, output);
    delete &alg;
}

void cb_dequantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.cb_dequantizate(tile, output);
    delete &alg;
}
void cr_dequantizate(const cv::Mat& tile, cv::Mat& output){
    AxDCT_algorithm& alg = *(new __USER_DEFAULT_ALGORITHM);
    alg.cr_dequantizate(tile, output);
    delete &alg;
}
// static const int __ZIGZAG[64] =
// {
//      0,  1,  8, 16,  9,  2,  3, 10,
//     17, 24, 32, 25, 18, 11,  4,  5,
//     12, 19, 26, 33, 40, 48, 41, 34,
//     27, 20, 13,  6,  7, 14, 21, 28,
//     35, 42, 49, 56, 57, 50, 43, 36,
//     29, 22, 15, 23, 30, 37, 44, 51,
//     58, 59, 52, 45, 38, 31, 39, 46,
//     53, 60, 61, 54, 47, 55, 62, 63,
// };

// void zigzag_encode(const cv::Mat& input, cv::Mat& output)
// {
//     cv::Mat temp(8,8,CV_64FC1);

//     for (int i=0; i<64; i++) buf [i] = data[__ZIGZAG[i]];
//     for (int i=0; i<64; i++) data[i] = buf[i];
// }

// void zigzag_decode(int *data)
// {
//     int buf[64], i;
//     for (i=0; i<64; i++) buf [__ZIGZAG[i]] = data[i];
//     for (i=0; i<64; i++) data[i] = buf[i];
// }