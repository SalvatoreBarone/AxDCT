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
 * @file   BC12.cpp
 * @author Andrea Aletto
 * @date   11 feb 2019
 * @brief  Implementation of BC12 algorithm class
 ******************************************************************************/

#include "BC12.h"

cv::Mat BC12::getT(){

    cv::Mat T = cv::Mat::zeros(8,8,CV_16S);

    T.at<int16_t>(0, 0) = 1; 
    T.at<int16_t>(0, 1) = 1;
    T.at<int16_t>(0, 2) = 1;
    T.at<int16_t>(0, 3) = 1;
    T.at<int16_t>(0, 4) = 1;
    T.at<int16_t>(0, 5) = 1;
    T.at<int16_t>(0, 6) = 1;
    T.at<int16_t>(0, 7) = 1;

    T.at<int16_t>(1, 0) = 1; 
    T.at<int16_t>(1, 7) = -1;

    T.at<int16_t>(2, 0) = 1;
    T.at<int16_t>(2, 3) = -1;
    T.at<int16_t>(2, 4) = -1;
    T.at<int16_t>(2, 7) = 1;

    T.at<int16_t>(3, 2) = -1;
    T.at<int16_t>(3, 5) = 1;

    T.at<int16_t>(4, 0) = 1; 
    T.at<int16_t>(4, 1) = -1;
    T.at<int16_t>(4, 2) = -1;
    T.at<int16_t>(4, 3) = 1;
    T.at<int16_t>(4, 4) = 1;
    T.at<int16_t>(4, 5) = -1;
    T.at<int16_t>(4, 6) = -1;
    T.at<int16_t>(4, 7) = 1;

    T.at<int16_t>(5, 1) = -1;
    T.at<int16_t>(5, 6) = 1;

    T.at<int16_t>(6, 1) = -1;
    T.at<int16_t>(6, 2) = 1;
    T.at<int16_t>(6, 5) = 1;
    T.at<int16_t>(6, 6) = -1;

    T.at<int16_t>(7, 3) = -1;
    T.at<int16_t>(7, 4) = 1;

    return T;
}

cv::Mat BC12::getD(){

    cv::Mat D = cv::Mat::zeros(8,8,CV_16S);

    D.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, -8> {1/sqrt(8)});
    D.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, -8> {1/sqrt(2)});
    D.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, -8> {0.5      });
    D.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, -8> {1/sqrt(2)});
    D.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, -8> {1/sqrt(8)});
    D.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, -8> {1/sqrt(2)});
    D.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, -8> {0.5      });
    D.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, -8> {1/sqrt(2)});

    return D;
}


/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ Q
 * 
 * This matrix should *multiply* the transformated luma image block 
 */ 
cv::Mat BC12::getLumaFullQuantMatrix(){

    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.007812500000000} );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.020833333333333} );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.012626906806903} );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017857142857143} );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006944444444444} );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.010416666666667} );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003607687659115} );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003472222222222} );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.022727272727273} );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.041666666666667} );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.027196414661021} );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.029411764705882} );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.011363636363636} );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.014285714285714} );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005524271728020} );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005434782608696} );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017677669529664} );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.025253813613805} );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.015625000000000} );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.016070608663331} );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004777748521531} );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006428243465332} );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003205128205128} );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003721614637824} );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.015625000000000} );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.026315789473684} );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.014731391274720} );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017241379310345} );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004464285714286} );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.007812500000000} );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004063832075785} );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005102040816327} );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005208333333333} );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.009615384615385} );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004419417382416} );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004901960784314} );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001838235294118} );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003086419753086} );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001716278595113} );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002232142857143} );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006250000000000} );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.008620689655172} );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006202691063040} );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005747126436782} );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002293577981651} );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004807692307692} );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002921928847878} );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005000000000000} );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003466209711699} );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005892556509888} );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003623188405797} );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004419417382416} );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001716278595113} );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003128791067197} );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002083333333333} );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003432557190226} );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004098360655738} );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.009090909090909} );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006313453403451} );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.008064516129032} );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003246753246753} );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005434782608696} );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003500528619735} );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005050505050505} );

    return luma_q_mat;
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ CQ
 * 
 * This matrix should *multiply* the transformated chroma image block 
 */ 
cv::Mat BC12::getChromaFullQuantMatrix(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013888888888889 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007365695637360 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013888888888889 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.023809523809524 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013598207330511 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007575757575758 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007365695637360 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013598207330511 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007575757575758 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );

    return chroma_q_mat;
}

cv::Mat BC12::getYQuantizationMatix(){
    return getLumaFullQuantMatrix();
}

cv::Mat BC12::getCbQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BC12::getCrQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BC12::getYDequantizationMatix(){
    return this->getStandardQ();
}

cv::Mat BC12::getCrDequantizationMatix(){
    return this->getStandardCQ();
}

cv::Mat BC12::getCbDequantizationMatix(){
    return this->getStandardCQ();
}

void BC12::dct1d(const cv::Mat& input, cv::Mat& output){
    /* This function implements the 1D transformation of the algorithm BC12:

       output = T * input  

       where:
            * T is the trannsformation matrix, of size 8x8
            * input is the input column vector to transform, ofize 8x1 
            * output is the transformated input vector, of size 8x1 
    */

    assert(( (input.rows == 8) && (input.cols==1) ) && "Column vector of size 8x1 is needed for 1D-DCT.");
    assert( (input.type() == CV_16S) && "Unable to compute AxDCT-1D: element of type CV_16S required.");

    /****** Transformation using matrix factorization X = P*A3*A2*A1 * x --- 14 adds + 4 sign inversions ******/
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
    int16_t x4b = x3a + (-x4a);
    int16_t x5b = x2a + (-x5a);
    int16_t x6b = x1a + (-x6a);
    int16_t x7b = x0a + (-x7a);

    int16_t x0c = x0b + x3b;
    int16_t x1c = x1b + x2b;
    int16_t x2c = x1b + (-x2b);
    int16_t x3c = x0b + (-x3b);
    int16_t x4c = -x4b;
    int16_t x5c = -x5b;
    int16_t x6c = -x6b;
    int16_t x7c = x7b;

    int16_t x0d = x0c + x1c;
    int16_t x1d = x0c + (-x1c);
    int16_t x2d = -x2c;
    int16_t x3d = x3c;
    int16_t x4d = x4c;
    int16_t x5d = x5c;
    int16_t x6d = x6c;
    int16_t x7d = x7c;

    output.at<int16_t>(0,0) = x0d;
    output.at<int16_t>(1,0) = x7d;
    output.at<int16_t>(2,0) = x3d;
    output.at<int16_t>(3,0) = x5d;
    output.at<int16_t>(4,0) = x1d;
    output.at<int16_t>(5,0) = x6d;
    output.at<int16_t>(6,0) = x2d;
    output.at<int16_t>(7,0) = x4d;
}