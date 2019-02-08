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
using namespace std;;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    Mat originalImage = imread( argv[1], IMREAD_GRAYSCALE );
    Mat image = imread( argv[1], IMREAD_GRAYSCALE );

    assert( originalImage.data && "No image data");

    namedWindow("Original Image", WINDOW_AUTOSIZE );
    imshow("Original Image", originalImage);

    //TODO: In image ho la matrice dei pixel dell'immagine

    // for(int i=0;i<512;i+=8){
    //     for(int j=0;j<512;j+=8){
    //         Mat block = image( cv::Rect(i, j, 8, 8) );
    //         block.at<uchar>(0,0) = 128;
    //     }
    // }

    //matrix_mult(originalImage, cv::Mat::zeros(512,512,CV_8U), image, 0);
    cerr<<__PRETTY_FUNCTION__<<__LINE__<<endl;
    Mat **v = splitInTiles(originalImage);
    cerr<<__PRETTY_FUNCTION__<<__LINE__<<endl;
    for(int i=0; i<2;i++){
        for(int j=0; j<2; j++){
            char name[255];
            sprintf(name, "Modified image (%d, %d)", i, j);
            namedWindow(name, WINDOW_AUTOSIZE );
            imshow(name, v[i][j]);
        }
    }

    // namedWindow("Modified Image", WINDOW_AUTOSIZE );
    // imshow("Modified Image", image);
    
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
                RES.at<unsigned char>(i,j) += A.at<unsigned char>(i,k) * B.at<unsigned char>(k,j);
            }
        }
    }    
}

cv::Mat **splitInTiles(const cv::Mat &input){
    if( ((input.rows%8) != 0) || ((input.cols%8) != 0) ) return NULL;
    int stride = 256;
    int r=input.rows/stride;
    int c=input.cols/stride;

    Mat **ret = new cv::Mat*[r];

    cerr<<__PRETTY_FUNCTION__<<__LINE__<<endl;
    for(int i=0; i<r; i++){
        ret[i] = new cv::Mat[c];
    }

    cerr<<__PRETTY_FUNCTION__<<__LINE__<<endl;

    for(int i=0; i<r; i++){
        for(int j=0; j<c; j++){
            (input(cv::Rect(i*stride,j*stride,stride,stride))).copyTo(ret[i][j]);
        }
    }

    return ret;
}