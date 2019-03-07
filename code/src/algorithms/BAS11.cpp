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
 * @file   BAS11.cpp
 * @author Andrea Aletto
 * @date   21 feb 2019
 * @brief  Implementation of BAS11 algorithm class
 ******************************************************************************/

#include "BAS11.h"
#include "user_defines.h"

static cv::Mat luma0Mat();
static cv::Mat luma05Mat();
static cv::Mat luma1Mat();
static cv::Mat luma2Mat();

static cv::Mat chroma0Mat();
static cv::Mat chroma05Mat();
static cv::Mat chroma1Mat();
static cv::Mat chroma2Mat();

BAS11::BAS11(double a) : AxDCT_algorithm() {
    assert(((a == 0.0) || (a == 1.0) || (a == 0.5) || (a == 2.0)) && "\n'a' param must be one of the following values: 0.0, 0.5, 1.0, 2.0");

    this->a = a;
}

cv::Mat BAS11::getD(){

    cv::Mat D = cv::Mat::zeros(8,8,CV_16S);

    D.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, -16> {   1/sqrt(8)              });
    D.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, -16> {   0.5                    });
    D.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, -16> {   1/sqrt( 4+4*(a*a) )    });
    D.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, -16> {   1/sqrt(2)              });
    D.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, -16> {   1/sqrt(8)              });
    D.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, -16> {   1/sqrt(2)              });
    D.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, -16> {   0.5                    });
    D.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, -16> {   1/sqrt( 4+4*(a*a) )    });

    return D;
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ Q
 * 
 * This matrix should *multiply* the transformated luma image block 
 */ 
cv::Mat BAS11::getLumaFullQuantMatrix(){

    if(this->a == 0.0) return luma0Mat();
    else if(this->a == 0.5) return luma05Mat();
    else if(this->a == 1.0) return luma1Mat();
    else if(this->a == 2.0) return luma2Mat();
    else return cv::Mat::zeros(8,8,CV_16S);
}

cv::Mat luma0Mat(){


    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.007812500000000} );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.014731391274720} );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.012626906806903} );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017857142857143} );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006944444444444} );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.010416666666667} );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003607687659115} );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002455231879120} );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.016070608663331} );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.020833333333333} );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.019230769230769} );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.020797258270193} );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.008035304331665} );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.010101525445522} );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003906250000000} );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002717391304348} );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017677669529664} );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017857142857143} );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.015625000000000} );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.016070608663331} );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004777748521531} );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006428243465332} );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003205128205128} );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002631578947368} );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.015625000000000} );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.018608073189120} );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.014731391274720} );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.017241379310345} );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004464285714286} );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.007812500000000} );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004063832075785} );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003607687659115} );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005208333333333} );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006799103665255} );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004419417382416} );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004901960784314} );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001838235294118} );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003086419753086} );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001716278595113} );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001578363350863} );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006250000000000} );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006095748113677} );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.006202691063040} );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005747126436782} );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002293577981651} );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004807692307692} );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002921928847878} );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003535533905933} );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003466209711699} );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004166666666667} );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003623188405797} );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004419417382416} );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.001716278595113} );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003128791067197} );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002083333333333} );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002427184466019} );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002897978611420} );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004545454545455} );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.004464285714286} );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.005702474041827} );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002295801237619} );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.003842971636883} );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002475247524752} );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{0.002525252525253} );

    return luma_q_mat;
   
}

cv::Mat luma05Mat(){


    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.014731391274720 } );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011293848786316 } );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017857142857143 } );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010416666666667 } );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003607687659115 } );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002196026152895 } );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.016070608663331 } );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020833333333333 } );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017200522903845 } );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020797258270193 } );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008035304331665 } );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010101525445522 } );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003906250000000 } );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002430508671195 } );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015811388300842 } );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015971914124998 } );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012500000000000 } );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.014373989364402 } );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004273348189417 } );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005749595745761 } );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002866753817307 } );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002105263157895 } );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015625000000000 } );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.018608073189120 } );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013176156917368 } );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017241379310345 } );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004063832075785 } );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003226813938947 } );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006799103665255 } );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003952847075210 } );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004901960784314 } );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001838235294118 } );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003086419753086 } );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001411731098289 } );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006250000000000 } );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006095748113677 } );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005547855544155 } );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005747126436782 } );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002293577981651 } );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004807692307692 } );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002921928847878 } );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003162277660168 } );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003466209711699 } );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004166666666667 } );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003240678228261 } );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004419417382416 } );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003128791067197 } );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002083333333333 } );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002170939783980 } );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002592030868990 } );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004065578140909 } );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571428571429 } );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005100447838981 } );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002053427052057 } );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003437258326270 } );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002213928690594 } );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );

    return luma_q_mat;
   
}

cv::Mat luma1Mat(){


    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.014731391274720 } );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008928571428571 } );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017857142857143 } );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010416666666667 } );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003607687659115 } );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001736111111111 } );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.016070608663331 } );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020833333333333 } );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.013598207330511 } );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020797258270193 } );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008035304331665 } );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010101525445522 } );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003906250000000 } );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001921485818442 } );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012500000000000 } );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.012626906806903 } );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011363636363636 } );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003378378378378 } );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004545454545455 } );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002266367888418 } );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001315789473684 } );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015625000000000 } );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.018608073189120 } );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010416666666667 } );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017241379310345 } );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004063832075785 } );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002551020408163 } );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006799103665255 } );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003125000000000 } );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004901960784314 } );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001838235294118 } );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003086419753086 } );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001116071428571 } );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006250000000000 } );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006095748113677 } );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004385964912281 } );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005747126436782 } );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002293577981651 } );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004807692307692 } );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002921928847878 } );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002500000000000 } );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003466209711699 } );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004166666666667 } );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002561981091256 } );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004419417382416 } );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003128791067197 } );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002083333333333 } );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002049180327869 } );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003214121732666 } );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002232142857143 } );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004032258064516 } );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001623376623377 } );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002717391304348 } );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001750264309868 } );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );

    return luma_q_mat;
   
}

cv::Mat luma2Mat(){


    cv::Mat luma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    luma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.014731391274720 } );
    luma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005646924393158 } );
    luma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017857142857143 } );
    luma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006944444444444 } );
    luma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010416666666667 } );
    luma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003607687659115 } );
    luma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001098013076447 } );
    luma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.016070608663331 } );
    luma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020833333333333 } );
    luma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008600261451922 } );
    luma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.020797258270193 } );
    luma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008035304331665 } );
    luma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.010101525445522 } );
    luma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003906250000000 } );
    luma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001215254335598 } );
    luma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007905694150421 } );
    luma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007985957062499 } );
    luma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003125000000000 } );
    luma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007186994682201 } );
    luma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002136674094708 } );
    luma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002874797872880 } );
    luma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001433376908654 } );
    luma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000526315789474 } );
    luma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.015625000000000 } );
    luma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.018608073189120 } );
    luma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006588078458684 } );
    luma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.017241379310345 } );
    luma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    luma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007812500000000 } );
    luma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004063832075785 } );
    luma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001613406969474 } );
    luma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    luma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006799103665255 } );
    luma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001976423537605 } );
    luma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004901960784314 } );
    luma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001838235294118 } );
    luma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003086419753086 } );
    luma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000705865549145 } );
    luma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006250000000000 } );
    luma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006095748113677 } );
    luma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002773927772078 } );
    luma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005747126436782 } );
    luma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002293577981651 } );
    luma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004807692307692 } );
    luma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002921928847878 } );
    luma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001581138830084 } );
    luma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003466209711699 } );
    luma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004166666666667 } );
    luma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001620339114130 } );
    luma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004419417382416 } );
    luma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001716278595113 } );
    luma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003128791067197 } );
    luma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002083333333333 } );
    luma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001085469891990 } );
    luma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001296015434495 } );
    luma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002032789070454 } );
    luma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000892857142857 } );
    luma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002550223919491 } );
    luma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001026713526029 } );
    luma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001718629163135 } );
    luma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001106964345297 } );
    luma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );

    return luma_q_mat;
   
}

/**
 * The following values are been computed with the formula:
 *      (diag(D) * diag(D)' ) ./ CQ
 * 
 * This matrix should *multiply* the transformated chroma image block 
 */ 
cv::Mat BAS11::getChromaFullQuantMatrix(){
    if(this->a == 0.0) return chroma0Mat();
    else if(this->a == 0.5) return chroma05Mat();
    else if(this->a == 1.0) return chroma1Mat();
    else if(this->a == 2.0) return chroma2Mat();
    else return cv::Mat::zeros(8,8,CV_16S);
}

cv::Mat chroma0Mat(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007365695637360 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011904761904762 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009615384615385 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007365695637360 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009615384615385 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004464285714286 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );

    return chroma_q_mat;
}

cv::Mat chroma05Mat(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006588078458684 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011904761904762 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008600261451922 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006588078458684 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.008600261451922 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571428571429 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003194219858756 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002258654522727 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002020202020202 } );

    return chroma_q_mat;
}

cv::Mat chroma1Mat(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011904761904762 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006799103665255 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005208333333333 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.006799103665255 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002232142857143 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );

    return chroma_q_mat;
}

cv::Mat chroma2Mat(){
    cv::Mat chroma_q_mat = cv::Mat::zeros(8,8,CV_16S);

    chroma_q_mat.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.007352941176471 } );
    chroma_q_mat.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003294039229342 } );
    chroma_q_mat.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.009820927516480 } );
    chroma_q_mat.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.011904761904762 } );
    chroma_q_mat.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004300130725961 } );
    chroma_q_mat.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003294039229342 } );
    chroma_q_mat.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.004300130725961 } );
    chroma_q_mat.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000892857142857 } );
    chroma_q_mat.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );
    chroma_q_mat.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005319148936170 } );
    chroma_q_mat.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005356869554444 } );
    chroma_q_mat.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001262626262626 } );
    chroma_q_mat.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.005050505050505 } );
    chroma_q_mat.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001785623184815 } );
    chroma_q_mat.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.003571246369629 } );
    chroma_q_mat.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.002525252525253 } );
    chroma_q_mat.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );
    chroma_q_mat.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000798554964689 } );
    chroma_q_mat.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001597109929378 } );
    chroma_q_mat.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.001129327261364 } );
    chroma_q_mat.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, _n_DECIMAL_>{ 0.000505050505051 } );

    return chroma_q_mat;
}

cv::Mat BAS11::getYQuantizationMatix(){
    return getLumaFullQuantMatrix();
}

cv::Mat BAS11::getCbQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BAS11::getCrQuantizationMatix(){
    return getChromaFullQuantMatrix();
}

cv::Mat BAS11::getYDequantizationMatix(){
    return this->getStandardQ();
}

cv::Mat BAS11::getCrDequantizationMatix(){
    return this->getStandardCQ();
}

cv::Mat BAS11::getCbDequantizationMatix(){
    return this->getStandardCQ();
}


void BAS11::dct1d(const cv::Mat& input, cv::Mat& output){

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
    int16_t x4c = x4b;
    int16_t x5c = x5b;
    int16_t x6c = x6b + x7b;
    int16_t x7c = x7b - x6b;

    int16_t x0d = x0c + x1c;
    int16_t x1d = x0c - x1c;
    int16_t x2d;
    int16_t x3d;
    int16_t x4d = x4c;
    int16_t x5d = x5c;
    int16_t x6d = x6c;
    int16_t x7d = x7c;

    // Operations for x2d and x3d are:
    //      x2d = a*x2c;
    //      x3d = a*x3c - x2c;

    if( a == 0) {
        x2d = 0;
        x3d = - x2c;

    } else if( a == 0.5 ){
        x2d = (x2c >> 1);
        x3d = (x3c >> 1) - x2c;

    } else if( a == 1 ){
        x2d = x2c;
        x3d = x3c - x2c;

    } else if( a == 2 ){
        x2d = (x2c << 1);
        x3d = (x3c << 1) - x2c;

    }

    output.at<int16_t>(0,0) = x0d;
    output.at<int16_t>(1,0) = x6d;
    output.at<int16_t>(2,0) = x2d;
    output.at<int16_t>(3,0) = x5d;
    output.at<int16_t>(4,0) = x1d;
    output.at<int16_t>(5,0) = x4d;
    output.at<int16_t>(6,0) = x7d;
    output.at<int16_t>(7,0) = x3d;

}
