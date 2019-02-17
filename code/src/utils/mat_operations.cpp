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
 * @file   mat_operations.cpp
 * @author Andrea Aletto
 * @date   16 feb 2019
 * @brief  Implementation of utils functions for image structures handling
 ******************************************************************************/

#include "mat_operations.h"

template<typename T>
void matrix_mult(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type){
    assert((A.cols == B.rows) && "Bad product multiplication");
    assert( ((type == CV_16S)||(type == CV_8U)||(type == CV_64FC1)) && "Type currently not supported" );
    assert(
        (( ( std::is_same<T, unsigned char>::value ) && ( type == CV_8U    ) ) ||
         ( ( std::is_same<T, short int    >::value ) && ( type == CV_16S   ) ) ||
         ( ( std::is_same<T, double       >::value ) && ( type == CV_64FC1 ) ) ) &&
         "Template type and requested destination type are incompatible"
    );

    /* Init matrix for calc */
    cv::Mat first(A.rows, A.cols, type);
    cv::Mat second(B.rows, B.cols, type);

    /* Type conversion if needed */
    if(A.type() != type) A.convertTo(first, type);
    else  A.copyTo(first);

    if(B.type() != type) B.convertTo(second, type);
    else B.copyTo(second);

    /* Init output matrix if it is NULL or of wrong size */
    if( !( RES.rows == A.rows && RES.cols == B.cols && RES.type() == type ) ) RES = cv::Mat::zeros(A.rows, B.cols, type);

    cv::Mat ret = cv::Mat::zeros(A.rows, B.cols, type);

    for(int i=0; i<A.rows; i++){
        for(int j=0; j<B.cols; j++){
            for(int k=0; k<A.cols; k++){
                // The operation is: RES[i][j] += A[i][k] * B[k][j];
                ret.at<T>(i,j) = ret.at<T>(i,j) + first.at<T>(i,k) * second.at<T>(k,j);
            }
        }
    }  
    ret.copyTo(RES);  
}

/* SPECIALIZATION TEMPLATE FUNCTIONS */
template void matrix_mult<unsigned char>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);
template void matrix_mult<int16_t>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);
template void matrix_mult<double>(cv::Mat const& A, cv::Mat const& B, cv::Mat& RES, int type);

cv::Mat **splitInTiles(const cv::Mat &input, int blockSize){
    if( ((input.rows%blockSize) != 0) || ((input.cols%blockSize) != 0) ) return NULL;

    int r=input.rows/blockSize;
    int c=input.cols/blockSize;

    cv::Mat **ret = new cv::Mat*[r];

    for(int i=0; i<r; i++){
        ret[i] = new cv::Mat[c];
    }

    for(int i=0; i<r; i++){
        for(int j=0; j<c; j++){
            (input(cv::Rect(i*blockSize,j*blockSize,blockSize,blockSize))).copyTo(ret[i][j]);
        }
    }

    return ret;
}

cv::Mat mergeTiles(cv::Mat **tiles, int imgWidth, int imgLength, int blockSize, bool deallocTiles){

    assert((imgWidth % blockSize == 0) && "Width must be multiple of block size");
    assert((imgLength % blockSize == 0) && "Length must be multiple of block size");

    cv::Mat ret(imgWidth, imgLength, CV_64FC1);

    for(int i=0; i<imgWidth/blockSize; i++){
        for(int j=0; j<imgLength/blockSize; j++){
            tiles[i][j].copyTo(ret(cv::Rect(i*blockSize,j*blockSize,blockSize,blockSize)));
            if(deallocTiles) tiles[i][j].deallocate();
        }
    }

    return ret;
}