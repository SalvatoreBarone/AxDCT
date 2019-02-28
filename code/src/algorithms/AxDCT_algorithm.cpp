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
 * @brief  AxDCT algorithm base class implementation
 ******************************************************************************/

#include "AxDCT_algorithm.h"

void AxDCT_algorithm::dct(const cv::Mat& tile, cv::Mat& output){

    /* Copy input tile ina temp variable in order to overwrite values */ 
    cv::Mat temp;
    tile.copyTo(temp);

    /* Row-wise transformation is required, but AxDCT-1D takes in input a column vector.
       So first transpose the tile, second apply AxDCT-1D with tile columns. 
     */ 
    transpose(temp, temp);

    /* Note: dct_col temporary variable is needed for build issues,
       the operation is: 
            temp.col(j) <- AxDCT-1D( temp.col(j) )    
     */
    for(int j=0; j<temp.cols; j++){
        cv::Mat dct_col = cv::Mat::zeros(8,1, CV_16S);
        this->dct1d(temp.col(j), dct_col);
        dct_col.copyTo(temp.col(j));
    }

    /* Now the column of temp are the AxDCT-1D of the rows of the original tile. */

    /* At this point a column-wise transformation is required, but the column of temp are
       the AxDCT-1D transformation of the rows. So a transposition is needed. 
     */
    transpose(temp, temp);

    /* Now it's possible to apply the AxDCT-1D transformation (identical to the former transformation). */
    for(int j=0; j<temp.cols; j++){
        cv::Mat dct_col = cv::Mat::zeros(8,1, CV_16S);
        this->dct1d(temp.col(j), dct_col);
        dct_col.copyTo(temp.col(j));
    }    

    /* In temp it's stored the AxDCT-2D of the input tile. Finally copy it to the output matrix */
    temp.copyTo(output);
}

void AxDCT_algorithm::y_quantizate(const cv::Mat& tile, cv::Mat& output){
    // std::cerr<<"\n\Calling y_quantizate\n";
    quantizate(tile, this->getYQuantizationMatix(), output);
}

void AxDCT_algorithm::cr_quantizate(const cv::Mat& tile, cv::Mat& output){
    quantizate(tile, this->getCrQuantizationMatix(), output);
}

void AxDCT_algorithm::cb_quantizate(const cv::Mat& tile, cv::Mat& output){
    quantizate(tile, this->getCbQuantizationMatix(), output);
}

void AxDCT_algorithm::quantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
    
    cv::Mat tileDCT;
    tile.convertTo(tileDCT, CV_16S);
    output = cv::Mat::zeros(tile.rows, tile.cols, CV_16S);
     
    for(int i=0; i<output.rows; i++){
        for(int j=0; j<output.cols; j++){
            int16_t q_val = (int16_t) (Q.at<int16_t>(i,j)); 
            int16_t x_val = (int16_t) (tileDCT.at<int16_t>(i,j));

            unsigned int prod = (x_val * q_val);
            if( prod & (1U << 17) ){
                prod = prod >> 16;
                prod++;
            } else {
                prod = prod >> 16;
            }
            
            output.at<int16_t>(i,j) = (int16_t)(prod );
        }
    }
}

void AxDCT_algorithm::y_dequantizate(const cv::Mat& tile, cv::Mat& output) {
    dequantizate(tile, this->getYDequantizationMatix(), output);
}

void AxDCT_algorithm::cr_dequantizate(const cv::Mat& tile, cv::Mat& output){
    dequantizate(tile, this->getCrDequantizationMatix(), output);
}

void AxDCT_algorithm::cb_dequantizate(const cv::Mat& tile, cv::Mat& output){
    dequantizate(tile, this->getCrDequantizationMatix(), output);
}

void AxDCT_algorithm::dequantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
    assert((tile.type() == CV_16S) && "Wrong tile type");
    assert((Q.type() == CV_16S) && "Wrong Q type");
    cv::Mat deq = tile.mul(Q);

    deq.copyTo(output);
}

cv::Mat AxDCT_algorithm::getStandardQ(){
    cv::Mat Q = cv::Mat::zeros(8,8,CV_16S);

    Q.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, 0>{16} );
    Q.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, 0>{12} );
    Q.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, 0>{14} );
    Q.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, 0>{14} );
    Q.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, 0>{18} );
    Q.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, 0>{24} );
    Q.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, 0>{49} );
    Q.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, 0>{72} );
    Q.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, 0>{11} );
    Q.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, 0>{12} );
    Q.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, 0>{13} );
    Q.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, 0>{17} );
    Q.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, 0>{22} );
    Q.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, 0>{35} );
    Q.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, 0>{64} );
    Q.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, 0>{92} );
    Q.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, 0>{10} );
    Q.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, 0>{14} );
    Q.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, 0>{16} );
    Q.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, 0>{22} );
    Q.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, 0>{37} );
    Q.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, 0>{55} );
    Q.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, 0>{78} );
    Q.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, 0>{95} );
    Q.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, 0>{16} );
    Q.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, 0>{19} );
    Q.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, 0>{24} );
    Q.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, 0>{29} );
    Q.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, 0>{56} );
    Q.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, 0>{64} );
    Q.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, 0>{87} );
    Q.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, 0>{98} );
    Q.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, 0>{24}  );
    Q.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, 0>{26}  );
    Q.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, 0>{40}  );
    Q.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, 0>{51}  );
    Q.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, 0>{68}  );
    Q.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, 0>{81}  );
    Q.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, 0>{103} );
    Q.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, 0>{112} );
    Q.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, 0>{40} );
    Q.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, 0>{58} );
    Q.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, 0>{57} );
    Q.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, 0>{87} );
    Q.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, 0>{109} );
    Q.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, 0>{104} );
    Q.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, 0>{121} );
    Q.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, 0>{100} );
    Q.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, 0>{51} );
    Q.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, 0>{60} );
    Q.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, 0>{69} );
    Q.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, 0>{80} );
    Q.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, 0>{103} );
    Q.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, 0>{113} );
    Q.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, 0>{120} );
    Q.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, 0>{103} );
    Q.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, 0>{61} );
    Q.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, 0>{55} );
    Q.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, 0>{56} );
    Q.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, 0>{62} );
    Q.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, 0>{77} );
    Q.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, 0>{92} );
    Q.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, 0>{101} );
    Q.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, 0>{99} );

    return Q;
}

cv::Mat AxDCT_algorithm::getStandardCQ(){
    cv::Mat CQ = cv::Mat::zeros(8,8,CV_16S);

    CQ.at<int16_t>(0, 0) = to_rep( fixed_point<int16_t, 0>{17});
    CQ.at<int16_t>(1, 0) = to_rep( fixed_point<int16_t, 0>{18});
    CQ.at<int16_t>(2, 0) = to_rep( fixed_point<int16_t, 0>{24});
    CQ.at<int16_t>(3, 0) = to_rep( fixed_point<int16_t, 0>{47});
    CQ.at<int16_t>(4, 0) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 0) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 0) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 0) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 1) = to_rep( fixed_point<int16_t, 0>{18});
    CQ.at<int16_t>(1, 1) = to_rep( fixed_point<int16_t, 0>{21});
    CQ.at<int16_t>(2, 1) = to_rep( fixed_point<int16_t, 0>{26});
    CQ.at<int16_t>(3, 1) = to_rep( fixed_point<int16_t, 0>{66});
    CQ.at<int16_t>(4, 1) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 1) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 1) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 1) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 2) = to_rep( fixed_point<int16_t, 0>{24});
    CQ.at<int16_t>(1, 2) = to_rep( fixed_point<int16_t, 0>{26});
    CQ.at<int16_t>(2, 2) = to_rep( fixed_point<int16_t, 0>{56});
    CQ.at<int16_t>(3, 2) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 2) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 2) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 2) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 2) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 3) = to_rep( fixed_point<int16_t, 0>{47});
    CQ.at<int16_t>(1, 3) = to_rep( fixed_point<int16_t, 0>{66});
    CQ.at<int16_t>(2, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(3, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 3) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(1, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(2, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(3, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 4) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(1, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(2, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(3, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 5) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(1, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(2, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(3, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 6) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(0, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(1, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(2, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(3, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(4, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(5, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(6, 7) = to_rep( fixed_point<int16_t, 0>{99});
    CQ.at<int16_t>(7, 7) = to_rep( fixed_point<int16_t, 0>{99});

    return CQ;
}
