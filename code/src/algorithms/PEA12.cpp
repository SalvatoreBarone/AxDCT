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

    cv::Mat D = cv::Mat::zeros(8,8,CV_16S);   

    D.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(2) });
    D.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(3) });
    D.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(5) });
    D.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(3) });
    D.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(2) });
    D.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(3) });
    D.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(5) });
    D.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, -16> { 2*sqrt(3) });

    return D;
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ Q
 * 
 * This matrix should *multiply* the transformated luma image block 
 */ 
cv::Mat PEA12::getLumaFullQuantMatrix(){

    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008505172717997 } );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005646924393158 } );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007290148043998 } );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004252586358999 } );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001613406969474 } );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001417528786333 } );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009278370237815 } );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004965363264368 } );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004901960784314 } );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004639185118908 } );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002380952380952 } );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001008589413075 } );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000905797101449 } );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007905694150421 } );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004610694459771 } );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003125000000000 } );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002934078292581 } );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002136674094708 } );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001173631317033 } );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000641025641026 } );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000679470762493 } );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006378879538498 } );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004385964912281 } );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002689571768200 } );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002873563218391 } );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001822537010999 } );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001302083333333 } );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000741950832607 } );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000850340136054 } );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003925464331383 } );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001976423537605 } );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002001217110117 } );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001838235294118 } );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001260025587851 } );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000767543121400 } );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000911268505500 } );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002551551815399 } );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001436781609195 } );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001132451270821 } );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000957854406130 } );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000936349290055 } );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000801282051282 } );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000533468780469 } );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000833333333333 } );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001550136107926 } );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001075828707280 } );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000724637681159 } );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000806871530460 } );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000767543121400 } );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000571236481742 } );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000416666666667 } );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000626696334338 } );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001673148731409 } );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001515151515152 } );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001152673614943 } );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001344086021505 } );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001325481462545 } );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000905797101449 } );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000639106162740 } );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );

    return luma_q_mat;
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ CQ
 * 
 * This matrix should *multiply* the transformated chroma image block 
 */ 
cv::Mat PEA12::getChromaFullQuantMatrix(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005670115145331 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003294039229342 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002171533459914 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005670115145331 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003968253968254 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002482681632184 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003294039229342 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002482681632184 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000892857142857 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002171533459914 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001030930026424 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000652017398351 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000841750841751 } );

    return chroma_q_mat;
}

cv::Mat PEA12::getYQuantizationMatix(){
    return getLumaFullQuantMatrix();
}

cv::Mat PEA12::getCbQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat PEA12::getCrQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat PEA12::getYDequantizationMatix(){
    return this->getStandardQ();
}

cv::Mat PEA12::getCrDequantizationMatix(){
    return this->getStandardCQ();
}

cv::Mat PEA12::getCbDequantizationMatix(){
    return this->getStandardCQ();
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
    int16_t x4c = x5b + x6b; 
    x4c = x4c + (x7b << 1);
    int16_t x5c = x7b - (x5b << 1); 
    x5c = x5c - x4b;
    int16_t x6c = x4b - (x6b << 1);
    x6c = x6c + x7b;
    int16_t x7c = x5b - (x4b << 1);
    x7c = x7c - x6b;

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
