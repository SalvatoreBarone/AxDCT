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

static void quantizate(const cv::Mat&, const cv::Mat&, cv::Mat&);
static void dequantizate(const cv::Mat&, const cv::Mat&, cv::Mat&);


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

    cv::Mat D = cv::Mat::zeros(8,8,CV_64FC1);

    D.at<double>(0, 0) = 1/sqrt(8);
    D.at<double>(1, 1) = 1/sqrt(2);
    D.at<double>(2, 2) = 0.5      ;
    D.at<double>(3, 3) = 1/sqrt(2);
    D.at<double>(4, 4) = 1/sqrt(8);
    D.at<double>(5, 5) = 1/sqrt(2);
    D.at<double>(6, 6) = 0.5      ;
    D.at<double>(7, 7) = 1/sqrt(2);   

    return D;
}

cv::Mat BC12::getQ(){

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

cv::Mat BC12::getCQ(){

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

cv::Mat BC12::getYQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getQ();
    return quantizationMatrix;

}

cv::Mat BC12::getCbQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getCQ();
    return quantizationMatrix;
}

cv::Mat BC12::getCrQuantizationMatix(){
    cv::Mat D_t, quantizationMatrix = cv::Mat::zeros(8,8, CV_64FC1);
    
    transpose(this->getD().diag(), D_t);
    matrix_mult<double>(this->getD().diag(), D_t, quantizationMatrix, CV_64FC1);

    quantizationMatrix /= this->getCQ();
    return quantizationMatrix;
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

    /***** DIRECT TRANSFORMATION X=T*x --- 24 adds *****/

    // int16_t x0 = input.at<int16_t>(0,0);
    // int16_t x1 = input.at<int16_t>(1,0);
    // int16_t x2 = input.at<int16_t>(2,0);
    // int16_t x3 = input.at<int16_t>(3,0);
    // int16_t x4 = input.at<int16_t>(4,0);
    // int16_t x5 = input.at<int16_t>(5,0);
    // int16_t x6 = input.at<int16_t>(6,0);
    // int16_t x7 = input.at<int16_t>(7,0);

    // output.at<int16_t>(0,0) = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;
    // output.at<int16_t>(1,0) = x0 + (-x7);
    // output.at<int16_t>(2,0) = x0 + (-x3) + (-x4) + x7;
    // output.at<int16_t>(3,0) = x5 + (-x2);
    // output.at<int16_t>(4,0) = x0 + (-x1) + (-x2) + x3 + x4 + (-x5) + (-x6) + x7;
    // output.at<int16_t>(5,0) = x6 + (-x1);
    // output.at<int16_t>(6,0) = x2 + (-x1) + x5 + (-x6);
    // output.at<int16_t>(7,0) = x4 + (-x3);



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

void BC12::y_quantizate(const cv::Mat& tile, cv::Mat& output){
    quantizate(tile, this->getYQuantizationMatix(), output);
}

void BC12::cr_quantizate(const cv::Mat& tile, cv::Mat& output){
    quantizate(tile, this->getCrQuantizationMatix(), output);
}

void BC12::cb_quantizate(const cv::Mat& tile, cv::Mat& output){
    quantizate(tile, this->getCbQuantizationMatix(), output);
}

void quantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
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

void BC12::y_dequantizate(const cv::Mat& tile, cv::Mat& output) {
    dequantizate(tile, this->getQ(), output);
}

void BC12::cr_dequantizate(const cv::Mat& tile, cv::Mat& output){
    dequantizate(tile, this->getCQ(), output);
}

void BC12::cb_dequantizate(const cv::Mat& tile, cv::Mat& output){
    dequantizate(tile, this->getCQ(), output);
}

void dequantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
    assert(tile.type() == CV_64FC1 && "Wrong tile type");
    assert(Q.type() == CV_64FC1 && "Wrong Q type");
    cv::Mat deq = tile.mul(Q);
    deq.copyTo(output);
}