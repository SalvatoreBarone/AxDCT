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

#define CHECKPOINT (std::cerr<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    cv::Mat bgrImg = imread( argv[1], cv::IMREAD_COLOR );
    cv::Mat ycrcbImg(512,512,CV_8U);

    assert( bgrImg.data && "No image data");

    cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE );
    imshow("Original Image", bgrImg);


    cv::cvtColor(bgrImg, ycrcbImg, cv::COLOR_BGR2YCrCb);
    std::vector<cv::Mat> chan(3);
    cv::split(ycrcbImg, chan);

    int blockSize = 8;
    cv::Mat **tiles = splitInTiles(chan[0], 8);

    cv::Mat T, D, Q;

    retrieveParameters(AxDCT_Algorithm::BC12, T, D, Q);
    
    cv::Mat res = cv::Mat::zeros(8,8,CV_16S);
    matrix_mult(T, cv::Mat::eye(8,8, CV_16S), res, CV_16S);
    PRINT_MAT(res, "T*I")

    for(int i=0;i<chan[0].rows/blockSize;i++){
        for(int j=0;j<chan[0].cols/blockSize;j++){
            if(i==0 && j==0) {
                PRINT_MAT(tiles[i][j], "original tile");
            }
            
            AxDCT(tiles[i][j], T, tiles[i][j]);
            if(i==0 && j==0) {
                PRINT_MAT(tiles[i][j], "DCT unquantizated");
            }

            quantizate(tiles[i][j], D, Q, tiles[i][j]);
            // tiles[i][j].convertTo(tiles[i][j], CV_8U);
            if(i==0 && j==0) {
                PRINT_MAT(tiles[i][j], "DCT quantizated");
            }

            if(i==0 && j==0) {
                PRINT_MAT(Q, "Luma Quantization Matrix");
                PRINT_MAT(T, "Transformation Matrix");
                PRINT_MAT(D, "Diagonal Matrix");
            }

            dequantizate(tiles[i][j], Q, tiles[i][j]);
            // tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
            if(i==0 && j==0) {
                PRINT_MAT(tiles[i][j], "DCT dequantizated");
            }

            cv::idct(tiles[i][j], tiles[i][j]);
            
            tiles[i][j].convertTo(tiles[i][j], CV_8U);
            if(i==0 && j==0) {
                PRINT_MAT(tiles[i][j], "IDCT");
            }
        }
    }
    
    
    
    (mergeTiles(tiles, chan[0].rows, chan[0].cols)).copyTo(chan[0]);
    chan[0].convertTo(chan[0], CV_8U);

    cv::merge(chan, ycrcbImg);
    cv::cvtColor(ycrcbImg, ycrcbImg, cv::COLOR_YCrCb2BGR);

    cv::namedWindow("Modified Image", cv::WINDOW_AUTOSIZE );
    imshow("Modified Image", ycrcbImg);
    
    cv::waitKey(0);
    return 0;
}

void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type){
    assert((A.cols == B.rows) && "Bad product multiplication");

    cv::Mat A_converted(A.rows, A.cols, type), B_converted(B.rows, B.cols, type);

    if(A.type() != type){
        A.convertTo(A_converted, type);
    } else {
        A.copyTo(A_converted);
    }

    if(B.type() != type){
        B.convertTo(B_converted, type);
    } else {
        B.copyTo(B_converted);
    }

    //TODO: check se res==NULL init RES
    cv::Mat ret = cv::Mat::zeros(A.rows, B.cols, type);

    for(int i=0; i<A.rows; i++){
        for(int j=0; j<B.cols; j++){
            for(int k=0; k<A.cols; k++){
                // RES[i][j] += A[i][k] * B[k][j];
                //TODO: parametrizzare il tipo
                if(type == CV_64FC1)
                    ret.at<double>(i,j) = ret.at<double>(i,j) + A_converted.at<double>(i,k) * B_converted.at<double>(k,j);
                else
                    ret.at<int16_t>(i,j) = ret.at<int16_t>(i,j) + A_converted.at<int16_t>(i,k) * B_converted.at<int16_t>(k,j);
            }
        }
    }  
    ret.copyTo(RES);  
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
    std::cout << ret.size();

    return ret;
}

void AxDCT(const cv::Mat& tile, const cv::Mat& T, cv::Mat& output){

    cv::Mat T_t;
    cv::transpose(T, T_t);
    matrix_mult(T, tile, output, CV_16S);
    matrix_mult(output, T_t, output, CV_16S);
}

void quantizate(const cv::Mat& tile, const cv::Mat& D, const cv::Mat& Q, cv::Mat& output){
    cv::Mat tileDCT;
    tile.convertTo(tileDCT, CV_64FC1);

    /*  D * (TXT) */
    matrix_mult(D, tileDCT, output, CV_64FC1);

    /*  (DTXT) * D' */
    matrix_mult(output, D, output, CV_64FC1);

    /*  (DTXTD)./Q */
    output /= Q;

    /* round */
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
    // tile.copyTo(output);
    // output *= Q;
    // output.mul(Q);
}


void retrieveParameters(const AxDCT_Algorithm alg, cv::Mat& T, cv::Mat& D, cv::Mat& Q){
    switch (alg)
    {
        case AxDCT_Algorithm::BC12 :
            T = cv::Mat::zeros(8,8,CV_16S);
            D = cv::Mat::zeros(8,8,CV_64FC1);
            Q = cv::Mat::zeros(8,8,CV_64FC1);

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
            

            D.at<double>(0, 0) = 1/sqrt(8);
            D.at<double>(1, 1) = 1/sqrt(2);
            D.at<double>(2, 2) = 0.5      ;
            D.at<double>(3, 3) = 1/sqrt(2);
            D.at<double>(4, 4) = 1/sqrt(8);
            D.at<double>(5, 5) = 1/sqrt(2);
            D.at<double>(6, 6) = 0.5      ;
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