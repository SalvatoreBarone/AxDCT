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
 * @file   BAS09.cpp
 * @author Andrea Aletto
 * @date   11 feb 2019
 * @brief  Implementation of BAS09 algorithm class
 ******************************************************************************/

#include "BAS09.h"

void BAS09::retrieveParameters(cv::Mat& T, cv::Mat& D, cv::Mat& Q, cv::Mat& CQ){
    
    T = BAS09::getT();
    D = BAS09::getD();
    Q = BAS09::getQ();
    CQ = BAS09::getCQ();
    
}

cv::Mat BAS09::getT(){

    cv::Mat T = cv::Mat::zeros(8,8,CV_16S);

    T.at<int16_t>(0, 0) = 1; 
    T.at<int16_t>(1, 0) = 1;
    T.at<int16_t>(2, 0) = 1;
    T.at<int16_t>(3, 0) = 1;
    T.at<int16_t>(4, 0) = 1;
    T.at<int16_t>(5, 0) = 1;
    T.at<int16_t>(6, 0) = 1;
    T.at<int16_t>(7, 0) = 1;

    T.at<int16_t>(0, 1) = 1; 
    T.at<int16_t>(1, 1) = 1;
    T.at<int16_t>(2, 1) = 1;
    T.at<int16_t>(3, 1) = -1;
    T.at<int16_t>(4, 1) = -1;
    T.at<int16_t>(5, 1) = -1;
    T.at<int16_t>(6, 1) = -1;
    T.at<int16_t>(7, 1) = -1;

    T.at<int16_t>(0, 2) = 1; 
    T.at<int16_t>(1, 2) = 1;
    T.at<int16_t>(2, 2) = -1;
    T.at<int16_t>(3, 2) = -1;
    T.at<int16_t>(4, 2) = -1;
    T.at<int16_t>(5, 2) = 1;
    T.at<int16_t>(6, 2) = 1;
    T.at<int16_t>(7, 2) = 1;

    T.at<int16_t>(0, 3) = 1; 
    T.at<int16_t>(1, 3) = 1;
    T.at<int16_t>(2, 3) = -1;
    T.at<int16_t>(3, 3) = -1;
    T.at<int16_t>(4, 3) = 1;
    T.at<int16_t>(5, 3) = 1;
    T.at<int16_t>(6, 3) = -1;
    T.at<int16_t>(7, 3) = -1;

    T.at<int16_t>(0, 4) = 1; 
    T.at<int16_t>(1, 4) = -1;
    T.at<int16_t>(2, 4) = -1;
    T.at<int16_t>(3, 4) = 1;
    T.at<int16_t>(4, 4) = 1;
    T.at<int16_t>(5, 4) = -1;
    T.at<int16_t>(6, 4) = -1;
    T.at<int16_t>(7, 4) = 1;

    T.at<int16_t>(0, 5) = 1; 
    T.at<int16_t>(1, 5) = -1;
    T.at<int16_t>(2, 5) = -1;
    T.at<int16_t>(3, 5) = 1;
    T.at<int16_t>(4, 5) = -1;
    T.at<int16_t>(5, 5) = -1;
    T.at<int16_t>(6, 5) = 1;
    T.at<int16_t>(7, 5) = -1;

    T.at<int16_t>(0, 6) = 1; 
    T.at<int16_t>(1, 6) = -1;
    T.at<int16_t>(2, 6) = 1;
    T.at<int16_t>(3, 6) = 1;
    T.at<int16_t>(4, 6) = -1;
    T.at<int16_t>(5, 6) = 1;
    T.at<int16_t>(6, 6) = -1;
    T.at<int16_t>(7, 6) = 1;

    T.at<int16_t>(0, 7) = 1; 
    T.at<int16_t>(1, 7) = -1;
    T.at<int16_t>(2, 7) = 1;
    T.at<int16_t>(3, 7) = -1;
    T.at<int16_t>(4, 7) = 1;
    T.at<int16_t>(5, 7) = -1;
    T.at<int16_t>(6, 7) = 1;
    T.at<int16_t>(7, 7) = -1;

    return T;
}

cv::Mat BAS09::getD(){

    cv::Mat D = cv::Mat::zeros(8,8,CV_64FC1);

    D.at<double>(0, 0) = 1/sqrt(8);
    D.at<double>(1, 1) = 0.5      ;
    D.at<double>(2, 2) = 1/sqrt(8);
    D.at<double>(3, 3) = 1/sqrt(2);
    D.at<double>(4, 4) = 1/sqrt(8);
    D.at<double>(5, 5) = 0.5      ;
    D.at<double>(6, 6) = 1/sqrt(8);
    D.at<double>(7, 7) = 1/sqrt(2);   

    return D;
}

cv::Mat BAS09::getQ(){

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

cv::Mat BAS09::getCQ(){

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