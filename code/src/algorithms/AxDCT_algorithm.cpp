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
    cv::Mat tileDCT;
    tile.convertTo(tileDCT, CV_64FC1);
    output = cv::Mat::zeros(tile.rows, tile.cols, CV_64FC1);
     
    for(int i=0; i<output.rows; i++){
        for(int j=0; j<output.cols; j++){
            output.at<double>(i,j) = round( tileDCT.at<double>(i,j) * Q.at<double>(i,j) );
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
    assert(tile.type() == CV_64FC1 && "Wrong tile type");
    assert(Q.type() == CV_64FC1 && "Wrong Q type");
    cv::Mat deq = tile.mul(Q);
    deq.copyTo(output);
}