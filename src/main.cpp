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

#define CHECKPOINT (cerr<<__PRETTY_FUNCTION__<<__LINE__<<endl)

using namespace cv;
using namespace std;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    Mat originalImage = imread( argv[1], IMREAD_GRAYSCALE );
    Mat image = imread( argv[1], IMREAD_GRAYSCALE );

    assert( originalImage.data && "No image data");

    namedWindow("Original Image", WINDOW_AUTOSIZE );
    imshow("Original Image", originalImage);

    //TODO: In image ho la matrice dei pixel dell'immagine

    //TODO: converto to ycbcr

    //matrix_mult(originalImage, cv::Mat::zeros(512,512,CV_8U), image, 0);
    int blockSize = 8;
    Mat **tiles = splitInTiles(originalImage, 8);

    Mat T;

    retrieveParameters(AxDCT_Algorithm::BC12, T);

    for(int i=0;i<originalImage.rows/blockSize;i++){
        for(int j=0;j<originalImage.cols/blockSize;j++){
            AxDCT(tiles[i][j], T, tiles[i][j]);
        }
    }
    
    for(int i=0;i<originalImage.rows/blockSize;i++){
        for(int j=0;j<originalImage.cols/blockSize;j++){
            image = mergeTiles(tiles, originalImage.rows, originalImage.cols);
        }
    }

    namedWindow("Modified Image", WINDOW_AUTOSIZE );
    imshow("Modified Image", image);
    
    waitKey(0);
    return 0;
}

void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type){
    assert((A.cols == B.rows) && "Bad product multiplication");

    RES = Mat(A.rows, B.cols, type);

    for(int i=0; i<A.rows; i++){
        for(int j=0; j<B.cols; j++){
            for(int k=0; k<A.cols; k++){
                // RES[i][j] += A[i][k] * B[k][j];
                RES.at<unsigned char>(i,j) = RES.at<unsigned char>(i,j) + A.at<unsigned char>(i,k) * B.at<unsigned char>(k,j);
            }
        }
    }    
}

cv::Mat **splitInTiles(const cv::Mat &input, int blockSize){
    if( ((input.rows%blockSize) != 0) || ((input.cols%blockSize) != 0) ) return NULL;

    int r=input.rows/blockSize;
    int c=input.cols/blockSize;

    Mat **ret = new cv::Mat*[r];

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

    cv::Mat ret(imgWidth, imgLength, CV_8U);

    for(int i=0; i<imgWidth/blockSize; i++){
        for(int j=0; j<imgLength/blockSize; j++){
            tiles[i][j].copyTo(ret(cv::Rect(i*blockSize,j*blockSize,blockSize,blockSize)));
            if(deallocTiles) tiles[i][j].deallocate();
        }
    }

    return ret;
}

void AxDCT(const cv::Mat& tile, const cv::Mat& T, cv::Mat& output){
    matrix_mult(T, tile, output);
    matrix_mult(output, T, output);
}

void retrieveParameters(const AxDCT_Algorithm alg, cv::Mat& T){
    switch (alg)
    {
        case AxDCT_Algorithm::BC12 :
            T = cv::Mat::zeros(8,8,CV_8U);

            T.at<unsigned char>(0, 0) = 1; 
            T.at<unsigned char>(0, 1) = 1;
            T.at<unsigned char>(0, 2) = 1;
            T.at<unsigned char>(0, 3) = 1;
            T.at<unsigned char>(0, 4) = 1;
            T.at<unsigned char>(0, 5) = 1;
            T.at<unsigned char>(0, 6) = 1;
            T.at<unsigned char>(0, 7) = 1;

            T.at<unsigned char>(1, 0) = 1; 
            T.at<unsigned char>(1, 7) = -1;

            T.at<unsigned char>(2, 0) = 1;
            T.at<unsigned char>(2, 3) = -1;
            T.at<unsigned char>(2, 4) = -1;
            T.at<unsigned char>(2, 7) = 1;

            T.at<unsigned char>(3, 2) = -1;
            T.at<unsigned char>(3, 5) = -1;

            T.at<unsigned char>(4, 0) = 1; 
            T.at<unsigned char>(4, 1) = -1;
            T.at<unsigned char>(4, 2) = -1;
            T.at<unsigned char>(4, 3) = 1;
            T.at<unsigned char>(4, 4) = 1;
            T.at<unsigned char>(4, 5) = -1;
            T.at<unsigned char>(4, 6) = -1;
            T.at<unsigned char>(4, 7) = 1;

            T.at<unsigned char>(5, 1) = -1;
            T.at<unsigned char>(5, 6) = 1;

            T.at<unsigned char>(6, 1) = -1;
            T.at<unsigned char>(6, 2) = 1;
            T.at<unsigned char>(6, 5) = 1;
            T.at<unsigned char>(6, 6) = -1;

            T.at<unsigned char>(7, 3) = -1;
            T.at<unsigned char>(7, 4) = 1;
            break;
    
        default:
            break;
    }
}