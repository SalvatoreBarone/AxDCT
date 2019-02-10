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
 * @file   main.cpp
 * @author Andrea Aletto
 * @date   4 feb 2019
 * @brief  Implementation of main executable functions
 ******************************************************************************/

#include "main.h"
#include <stdio.h>

#define CHECKPOINT (std::cerr<<__PRETTY_FUNCTION__<<__LINE__<<std::endl)
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    // Load img
    cv::Mat bgrImg = imread( argv[1], cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");

    // Declare an empty for dst image
    cv::Mat ycrcbImg;

    // Convert to ycrcb
    cv::cvtColor(bgrImg, ycrcbImg, cv::COLOR_BGR2YCrCb);

    // Split bgr into 3 channels
    cv::Mat bgrChan[3];
    cv::split(bgrImg, bgrChan);

    // Split ycrcb into 3 channels
    cv::Mat ycrcbChan[3];
    cv::split(ycrcbImg, ycrcbChan);

    // Print first block for each channel
    PRINT_MAT(ycrcbChan[0](cv::Rect(0, 0, 8, 8)), "LUMA (first 8x8 block)")
    PRINT_MAT(ycrcbChan[1](cv::Rect(0, 0, 8, 8)), "Cr (first 8x8 block)")
    PRINT_MAT(ycrcbChan[2](cv::Rect(0, 0, 8, 8)), "Cb (first 8x8 block)")

    PRINT_MAT(bgrChan[0](cv::Rect(0, 0, 8, 8)), "Blue (first 8x8 block)")
    PRINT_MAT(bgrChan[1](cv::Rect(0, 0, 8, 8)), "Green (first 8x8 block)")
    PRINT_MAT(bgrChan[2](cv::Rect(0, 0, 8, 8)), "Red (first 8x8 block)")

    // Show results 
    cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE );
    imshow("Original Image", bgrImg);

    cv::namedWindow("Converted Image", cv::WINDOW_AUTOSIZE );
    imshow("Converted Image", ycrcbImg);
    
    cv::waitKey(0);
    return 0;
}

void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type){
    assert((A.cols == B.rows) && "Bad product multiplication");

    RES = cv::Mat(A.rows, B.cols, type);

    for(int i=0; i<A.rows; i++){
        for(int j=0; j<B.cols; j++){
            for(int k=0; k<A.cols; k++){
                // RES[i][j] += A[i][k] * B[k][j];
                RES.at<double>(i,j) = RES.at<double>(i,j) + A.at<double>(i,k) * B.at<double>(k,j);
            }
        }
    }    
}

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

void AxDCT(const cv::Mat& tile, const cv::Mat& T, cv::Mat& output){

    cv::Mat T_t;
    cv::transpose(T, T_t);
    matrix_mult(T, tile-128, output);
    matrix_mult(output, T_t, output);
}

void quantizate(const cv::Mat& tile, const cv::Mat& D, const cv::Mat& Q, cv::Mat& output){
    cv::Mat tileDCT;
    tile.convertTo(tileDCT, CV_64FC1);
    matrix_mult(D, tileDCT, output, CV_64FC1);
    matrix_mult(output, D, output, CV_64FC1);
    output.mul(1/Q);
}

void dequantizate(const cv::Mat& tile, const cv::Mat& Q, cv::Mat& output){
    tile.copyTo(output);
    output.mul(Q);
}


void retrieveParameters(const AxDCT_Algorithm alg, cv::Mat& T, cv::Mat& D, cv::Mat& Q){
    switch (alg)
    {
        case AxDCT_Algorithm::BC12 :
            T = cv::Mat::zeros(8,8,CV_64FC1);
            D = cv::Mat::zeros(8,8,CV_64FC1);
            Q = cv::Mat::zeros(8,8,CV_64FC1);

            T.at<double>(0, 0) = 1; 
            T.at<double>(0, 1) = 1;
            T.at<double>(0, 2) = 1;
            T.at<double>(0, 3) = 1;
            T.at<double>(0, 4) = 1;
            T.at<double>(0, 5) = 1;
            T.at<double>(0, 6) = 1;
            T.at<double>(0, 7) = 1;

            T.at<double>(1, 0) = 1; 
            T.at<double>(1, 7) = -1;

            T.at<double>(2, 0) = 1;
            T.at<double>(2, 3) = -1;
            T.at<double>(2, 4) = -1;
            T.at<double>(2, 7) = 1;

            T.at<double>(3, 2) = -1;
            T.at<double>(3, 5) = 1;

            T.at<double>(4, 0) = 1; 
            T.at<double>(4, 1) = -1;
            T.at<double>(4, 2) = -1;
            T.at<double>(4, 3) = 1;
            T.at<double>(4, 4) = 1;
            T.at<double>(4, 5) = -1;
            T.at<double>(4, 6) = -1;
            T.at<double>(4, 7) = 1;

            T.at<double>(5, 1) = -1;
            T.at<double>(5, 6) = 1;

            T.at<double>(6, 1) = -1;
            T.at<double>(6, 2) = 1;
            T.at<double>(6, 5) = 1;
            T.at<double>(6, 6) = -1;

            T.at<double>(7, 3) = -1;
            T.at<double>(7, 4) = 1;
            

            D.at<double>(0, 0) = 1/sqrt(8);
            D.at<double>(1, 1) = 1/sqrt(2);
            D.at<double>(2, 2) = 1/2      ;
            D.at<double>(3, 3) = 1/sqrt(2);
            D.at<double>(4, 4) = 1/sqrt(8);
            D.at<double>(5, 5) = 1/sqrt(2);
            D.at<double>(6, 6) = 1/2      ;
            D.at<double>(7, 7) = 1/sqrt(2);    

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

            break;
    
        default:
            break;
    }
}