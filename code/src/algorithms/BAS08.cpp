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
 * @file   BAS08.cpp
 * @author Andrea Aletto
 * @date   11 feb 2019
 * @brief  Implementation of BAS08 algorithm class
 ******************************************************************************/

#include "BAS08.h"

cv::Mat BAS08::getD(){

    cv::Mat D = cv::Mat::zeros(8,8,CV_64FC1);

    D.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(8) });
    D.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(2) });
    D.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(5) });
    D.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(2) });
    D.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(8) });
    D.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, -8> { 0.5       });
    D.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(5) });
    D.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, -8> { 1/sqrt(2) });

    return D;
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ Q
 * 
 * This matrix should *multiply* the transformated luma image block 
 */ 
cv::Mat BAS08::getLumaFullQuantMatrix(){

    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020833333333333 } );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011293848786316 } );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017857142857143 } );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007365695637360 } );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003226813938947 } );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003472222222222 } );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.022727272727273 } );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.041666666666667 } );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.024325212770526 } );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.029411764705882 } );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011363636363636 } );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010101525445522 } );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004941058844013 } );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005434782608696 } );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015811388300842 } );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.022587697572631 } );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012500000000000 } );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.014373989364402 } );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004273348189417 } );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004065578140909 } );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002564102564103 } );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003328713326493 } );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015625000000000 } );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.026315789473684 } );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013176156917368 } );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017241379310345 } );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005524271728020 } );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003634801908240 } );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005102040816327 } );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009615384615385 } );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003952847075210 } );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004901960784314 } );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001838235294118 } );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002182428336996 } );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001535086242800 } );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002232142857143 } );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004419417382416 } );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006095748113677 } );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003922926276315 } );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004063832075785 } );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001621804544006 } );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002403846153846 } );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001847990064049 } );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003535533905933 } );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003100272215851 } );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005270462766947 } );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002898550724638 } );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003952847075210 } );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001535086242800 } );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001978821219026 } );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001666666666667 } );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003070172485600 } );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004098360655738 } );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009090909090909 } );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005646924393158 } );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008064516129032 } );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003246753246753 } );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003842971636883 } );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003130967980365 } );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );

    return luma_q_mat;

}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ CQ
 * 
 * This matrix should *multiply* the transformated chroma image block 
 */ 
cv::Mat BAS08::getChromaFullQuantMatrix(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013888888888889 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006588078458684 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013888888888889 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.023809523809524 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012162606385263 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007575757575758 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006588078458684 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012162606385263 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571428571429 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007575757575758 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );

    return chroma_q_mat;
}

cv::Mat BAS08::getYQuantizationMatix(){
    return getLumaFullQuantMatrix();
}

cv::Mat BAS08::getCbQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BAS08::getCrQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BAS08::getYDequantizationMatix(){
    return this->getStandardQ();
}

cv::Mat BAS08::getCrDequantizationMatix(){
    return this->getStandardCQ();
}

cv::Mat BAS08::getCbDequantizationMatix(){
    return this->getStandardCQ();
}


void BAS08::dct1d(const cv::Mat& input, cv::Mat& output){

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
    int16_t x1c = x6b + x7b;
    int16_t x2c = x1b + x2b;
    int16_t x3c = -x5b;
    int16_t x4c = x1b - x2b;
    int16_t x5c = x7b - x6b;
    int16_t x6c = x0b - x3b;
    int16_t x7c = -x4b;

    output.at<int16_t>(0,0) = x0c + x2c;
    output.at<int16_t>(1,0) = x1c;
    output.at<int16_t>(2,0) = (x4c >> 1) + x6c;
    output.at<int16_t>(3,0) = x3c;
    output.at<int16_t>(4,0) = x0c - x2c;
    output.at<int16_t>(5,0) = x5c;
    output.at<int16_t>(6,0) = (x6c >> 1) - x4c;
    output.at<int16_t>(7,0) = x7c;

}
