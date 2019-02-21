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
 * @file   PEA12.cpp
 * @author Andrea Aletto
 * @date   21 feb 2019
 * @brief  Implementation of PEA12 algorithm class
 ******************************************************************************/

#include "PEA12.h"

cv::Mat PEA12::getD(){

    cv::Mat D = cv::Mat::zeros(8,8,CV_64FC1);

    D.at<double>(0, 0) = 1/(2*sqrt(2));
    D.at<double>(1, 1) = 1/(2*sqrt(3));
    D.at<double>(2, 2) = 1/(2*sqrt(5));
    D.at<double>(3, 3) = 1/(2*sqrt(3));
    D.at<double>(4, 4) = 1/(2*sqrt(2));
    D.at<double>(5, 5) = 1/(2*sqrt(3));
    D.at<double>(6, 6) = 1/(2*sqrt(5));
    D.at<double>(7, 7) = 1/(2*sqrt(3));   

    return D;
}

cv::Mat PEA12::getQ(){

    cv::Mat Q = cv::Mat::zeros(8,8,CV_64FC1);

    Q.at<double>(0, 0) = 16;
    Q.at<double>(1, 0) = 12;
    Q.at<double>(2, 0) = 14;
    Q.at<double>(3, 0) = 14;
    Q.at<double>(4, 0) = 18;
    Q.at<double>(5, 0) = 24;
    Q.at<double>(6, 0) = 49;
    Q.at<double>(7, 0) = 72;

    Q.at<double>(0, 1) = 11;
    Q.at<double>(1, 1) = 12;
    Q.at<double>(2, 1) = 13;
    Q.at<double>(3, 1) = 17;
    Q.at<double>(4, 1) = 22;
    Q.at<double>(5, 1) = 35;
    Q.at<double>(6, 1) = 64;
    Q.at<double>(7, 1) = 92;

    Q.at<double>(0, 2) = 10;
    Q.at<double>(1, 2) = 14;
    Q.at<double>(2, 2) = 16;
    Q.at<double>(3, 2) = 22;
    Q.at<double>(4, 2) = 37;
    Q.at<double>(5, 2) = 55;
    Q.at<double>(6, 2) = 78;
    Q.at<double>(7, 2) = 95;

    Q.at<double>(0, 3) = 16;
    Q.at<double>(1, 3) = 19;
    Q.at<double>(2, 3) = 24;
    Q.at<double>(3, 3) = 29;
    Q.at<double>(4, 3) = 56;
    Q.at<double>(5, 3) = 64;
    Q.at<double>(6, 3) = 87;
    Q.at<double>(7, 3) = 98;

    Q.at<double>(0, 4) = 24;
    Q.at<double>(1, 4) = 26;
    Q.at<double>(2, 4) = 40;
    Q.at<double>(3, 4) = 51;
    Q.at<double>(4, 4) = 68;
    Q.at<double>(5, 4) = 81;
    Q.at<double>(6, 4) = 103;
    Q.at<double>(7, 4) = 112;

    Q.at<double>(0, 5) = 40;
    Q.at<double>(1, 5) = 58;
    Q.at<double>(2, 5) = 57;
    Q.at<double>(3, 5) = 87;
    Q.at<double>(4, 5) = 109;
    Q.at<double>(5, 5) = 104;
    Q.at<double>(6, 5) = 121;
    Q.at<double>(7, 5) = 100;

    Q.at<double>(0, 6) = 51;
    Q.at<double>(1, 6) = 60;
    Q.at<double>(2, 6) = 69;
    Q.at<double>(3, 6) = 80;
    Q.at<double>(4, 6) = 103;
    Q.at<double>(5, 6) = 113;
    Q.at<double>(6, 6) = 120;
    Q.at<double>(7, 6) = 103;

    Q.at<double>(0, 7) = 61;
    Q.at<double>(1, 7) = 55;
    Q.at<double>(2, 7) = 56;
    Q.at<double>(3, 7) = 62;
    Q.at<double>(4, 7) = 77;
    Q.at<double>(5, 7) = 92;
    Q.at<double>(6, 7) = 101;
    Q.at<double>(7, 7) = 99;

    return Q;
}

cv::Mat PEA12::getCQ(){

    cv::Mat CQ = cv::Mat::zeros(8,8,CV_64FC1);

    CQ.at<double>(0, 0) = 17;
    CQ.at<double>(1, 0) = 18;
    CQ.at<double>(2, 0) = 24;
    CQ.at<double>(3, 0) = 47;
    CQ.at<double>(4, 0) = 99;
    CQ.at<double>(5, 0) = 99;
    CQ.at<double>(6, 0) = 99;
    CQ.at<double>(7, 0) = 99;

    CQ.at<double>(0, 1) = 18;
    CQ.at<double>(1, 1) = 21;
    CQ.at<double>(2, 1) = 26;
    CQ.at<double>(3, 1) = 66;
    CQ.at<double>(4, 1) = 99;
    CQ.at<double>(5, 1) = 99;
    CQ.at<double>(6, 1) = 99;
    CQ.at<double>(7, 1) = 99;

    CQ.at<double>(0, 2) = 24;
    CQ.at<double>(1, 2) = 26;
    CQ.at<double>(2, 2) = 56;
    CQ.at<double>(3, 2) = 99;
    CQ.at<double>(4, 2) = 99;
    CQ.at<double>(5, 2) = 99;
    CQ.at<double>(6, 2) = 99;
    CQ.at<double>(7, 2) = 99;

    CQ.at<double>(0, 3) = 47;
    CQ.at<double>(1, 3) = 66;
    CQ.at<double>(2, 3) = 99;
    CQ.at<double>(3, 3) = 99;
    CQ.at<double>(4, 3) = 99;
    CQ.at<double>(5, 3) = 99;
    CQ.at<double>(6, 3) = 99;
    CQ.at<double>(7, 3) = 99;

    CQ.at<double>(0, 4) = 99;
    CQ.at<double>(1, 4) = 99;
    CQ.at<double>(2, 4) = 99;
    CQ.at<double>(3, 4) = 99;
    CQ.at<double>(4, 4) = 99;
    CQ.at<double>(5, 4) = 99;
    CQ.at<double>(6, 4) = 99;
    CQ.at<double>(7, 4) = 99;

    CQ.at<double>(0, 5) = 99;
    CQ.at<double>(1, 5) = 99;
    CQ.at<double>(2, 5) = 99;
    CQ.at<double>(3, 5) = 99;
    CQ.at<double>(4, 5) = 99;
    CQ.at<double>(5, 5) = 99;
    CQ.at<double>(6, 5) = 99;
    CQ.at<double>(7, 5) = 99;

    CQ.at<double>(0, 6) = 99;
    CQ.at<double>(1, 6) = 99;
    CQ.at<double>(2, 6) = 99;
    CQ.at<double>(3, 6) = 99;
    CQ.at<double>(4, 6) = 99;
    CQ.at<double>(5, 6) = 99;
    CQ.at<double>(6, 6) = 99;
    CQ.at<double>(7, 6) = 99;

    CQ.at<double>(0, 7) = 99;
    CQ.at<double>(1, 7) = 99;
    CQ.at<double>(2, 7) = 99;
    CQ.at<double>(3, 7) = 99;
    CQ.at<double>(4, 7) = 99;
    CQ.at<double>(5, 7) = 99;
    CQ.at<double>(6, 7) = 99;
    CQ.at<double>(7, 7) = 99;

    return CQ;
}

cv::Mat PEA12::getYQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getQ();
    return quantizationMatrix;

}

cv::Mat PEA12::getCbQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getCQ();
    return quantizationMatrix;
}

cv::Mat PEA12::getCrQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getCQ();
    return quantizationMatrix;
}

cv::Mat PEA12::getYDequantizationMatix(){
    return this->getQ();
}

cv::Mat PEA12::getCrDequantizationMatix(){
    return this->getCQ();
}

cv::Mat PEA12::getCbDequantizationMatix(){
    return this->getCQ();
}

void PEA12::dct1d(const cv::Mat& input, cv::Mat& output){

    assert(( (input.rows == 8) && (input.cols==1) ) && "Column vector of size 8x1 is needed for 1D-DCT.");
    assert( (input.type() == CV_16S) && "Unable to compute AxDCT-1D: element of type CV_16S required.");

    int16_t x0a = input.at<int16_t>(0,0);
    int16_t x1a = input.at<int16_t>(1,0);
    int16_t x2a = input.at<int16_t>(2,0);
    int16_t x3a = input.at<int16_t>(3,0);
    int16_t x4a = input.at<int16_t>(4,0);
    int16_t x5a = input.at<int16_t>(5,0);
    int16_t x6a = input.at<int16_t>(6,0);
    int16_t x7a = input.at<int16_t>(7,0);

    int16_t x0b = x0a + x7a;
    int16_t x1b = x1a + x6a;
    int16_t x2b = x2a + x5a;
    int16_t x3b = x3a + x4a;
    int16_t x4b = x3a - x4a;
    int16_t x5b = x2a - x5a;
    int16_t x6b = x1a - x6a;
    int16_t x7b = x0a - x7a;

    int16_t x0c = x0b + x3b;
    int16_t x1c = x1b + x2b;
    int16_t x2c = x1b - x2b;
    int16_t x3c = x0b - x3b;
    int16_t x4c = x5b + x6b + (x7b << 1);
    int16_t x5c = x7b - (x5b << 1) - x4b;
    int16_t x6c = x4b - (x6b << 1) + x7b;
    int16_t x7c = x5b - (x4b << 1) - x6b;

    int16_t x0d = x0c + x1c;
    int16_t x1d = x0c - x1c;
    int16_t x2d = x2c + (x3c << 1);
    int16_t x3d = x3c - (x2c << 1);
    int16_t x4d = x4c;
    int16_t x5d = x5c;
    int16_t x6d = x6c;
    int16_t x7d = x7c;

    output.at<int16_t>(0,0) = x0d;
    output.at<int16_t>(1,0) = x4d;
    output.at<int16_t>(2,0) = x2d;
    output.at<int16_t>(3,0) = x5d;
    output.at<int16_t>(4,0) = x1d;
    output.at<int16_t>(5,0) = x6d;
    output.at<int16_t>(6,0) = x3d;
    output.at<int16_t>(7,0) = x7d;

}
