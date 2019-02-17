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

static void AxDCT1D(const cv::Mat&, cv::Mat&);

void AxDCT(const cv::Mat& tile, cv::Mat& output){

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
        AxDCT1D(temp.col(j), dct_col);
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
        AxDCT1D(temp.col(j), dct_col);
        dct_col.copyTo(temp.col(j));
    }    

    /* In temp it's stored the AxDCT-2D of the input tile. Finally copy it to the output matrix */
    temp.copyTo(output);

}

void AxDCT1D(const cv::Mat& input, cv::Mat& output){
    /* This function implements the 1D transformation of the algorithm BC12:

       output = T * input  

       where:
            * T is the trannsformation matrix, of size 8x8
            * input is the input column vector to transform, ofize 8x1 
            * output is the transformated input vector, of size 8x1 
    */

    assert(( (input.rows == 8) && (input.cols==1) ) && "Column vector of size 8x1 is needed for 1D-DCT.");
    assert( (input.type() == CV_16S) && "Unable to compute AxDCT-1D: element of type CV_16S required.");

    /***** DIRECT TRANSFORMATION X=T*x --- 24 adds *****/
/*
    int16_t x0 = input.at<int16_t>(0,0);
    int16_t x1 = input.at<int16_t>(1,0);
    int16_t x2 = input.at<int16_t>(2,0);
    int16_t x3 = input.at<int16_t>(3,0);
    int16_t x4 = input.at<int16_t>(4,0);
    int16_t x5 = input.at<int16_t>(5,0);
    int16_t x6 = input.at<int16_t>(6,0);
    int16_t x7 = input.at<int16_t>(7,0);

    output.at<int16_t>(0,0) = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;
    output.at<int16_t>(1,0) = x0 + (-x7);
    output.at<int16_t>(2,0) = x0 + (-x3) + (-x4) + x7;
    output.at<int16_t>(3,0) = x5 + (-x2);
    output.at<int16_t>(4,0) = x0 + (-x1) + (-x2) + x3 + x4 + (-x5) + (-x6) + x7;
    output.at<int16_t>(5,0) = x6 + (-x1);
    output.at<int16_t>(6,0) = x2 + (-x1) + x5 + (-x6);
    output.at<int16_t>(7,0) = x4 + (-x3);
*/


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

void quantizate(const cv::Mat& tile, const cv::Mat& D, const cv::Mat& Q, cv::Mat& output){
    /*
        D matrix is merged with the Q matrix. Since D is diagonal, then: 
        
            D * (tile) * D' = (diag(D) * diag(D)') .* (tile)

        This result should be divided elem-wise by Q.

        An equivalent quantization is the following:
        F = tile .* Y
        where F is the DCT quantizated transformed tile and Y is the following matrix:

        Y = ( diag(D) * diag(D)' ) ./ Q

        Note that the Y matrix can be computed only once, offline.
        
    */
    cv::Mat tileDCT, D_t, DDt;
    tile.convertTo(tileDCT, CV_64FC1);
    
    transpose(D.diag(), D_t);
    matrix_mult<double>(D.diag(), D_t, DDt, CV_64FC1);
    // matrix_mult(D_diag, D_t, DDt, CV_64FC1);
    D_t.deallocate();   //cleanup

    output = tileDCT.mul(DDt);
    DDt.deallocate(); //cleanup
    output /= Q;

    /* round */ //TODO: check if this can be done in a non pixel-by-pixel way
    for(int i=0; i<output.rows; i++){
        for(int j=0; j<output.cols; j++){
            output.at<double>(i,j) = round(output.at<double>(i,j));
        }
    }
}

void dequantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
    assert(tile.type() == CV_64FC1 && "Wrong tile type");
    assert(Q.type() == CV_64FC1 && "Wrong Q type");
    cv::Mat deq = tile.mul(Q);
    deq.copyTo(output);
}