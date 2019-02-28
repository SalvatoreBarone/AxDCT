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

void AxDCT(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->dct(tile, output);
}

void y_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->y_quantizate(tile, output);
}

void cr_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->cr_quantizate(tile, output);
}

void cb_quantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->cb_quantizate(tile, output);
}

void y_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->y_dequantizate(tile, output);
}

void cb_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->cb_dequantizate(tile, output);
}

void cr_dequantizate(const cv::Mat& tile, cv::Mat& output, AxDCT_algorithm *alg){
    alg->cr_dequantizate(tile, output);
}

void transformImage(const cv::Mat& bgrImg, cv::Mat& output, AxDCT_algorithm *alg){

    // Declare an empty for dst image
    cv::Mat ycrcbImg;

    /* Convert BGR image into YCrCb */
    cv::cvtColor(bgrImg, ycrcbImg, cv::COLOR_BGR2YCrCb);

    /* Split YCbCr into 3 channels */
    std::vector<cv::Mat> chan(3);
    cv::split(ycrcbImg, chan);

    chan[0].convertTo(chan[0], CV_16S);
    chan[1].convertTo(chan[1], CV_16S);
    chan[2].convertTo(chan[2], CV_16S);

    /* Parameters for transformation */
    const int blockSize = 8;

    /********* LUMA *********/

    /* Split channel in blocks 8x8 */
    cv::Mat **tiles = splitInTiles(chan[0], blockSize);

    for(int i=0;i<chan[0].rows/blockSize;i++){
        for(int j=0;j<chan[0].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], tiles[i][j], alg);

            /* Quantization step */
            y_quantizate(tiles[i][j], tiles[i][j], alg);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[0].rows, chan[0].cols)).copyTo(chan[0]);

    /**********************/

    /********* Cr *********/

    /* Split channel in blocks 8x8 */
    tiles = splitInTiles(chan[1], blockSize);

    for(int i=0;i<chan[1].rows/blockSize;i++){
        for(int j=0;j<chan[1].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], tiles[i][j], alg);

            /* Quantization step */
            cr_quantizate(tiles[i][j], tiles[i][j], alg);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[1].rows, chan[1].cols)).copyTo(chan[1]);

    /**********************/

    /********* Cb *********/

    /* Split channel in blocks 8x8 */
    tiles = splitInTiles(chan[2], blockSize);

    for(int i=0;i<chan[2].rows/blockSize;i++){
        for(int j=0;j<chan[2].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], tiles[i][j], alg);

            /* Quantization step */
            cb_quantizate(tiles[i][j], tiles[i][j], alg);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[2].rows, chan[2].cols)).copyTo(chan[2]);

    /**********************/
    
    /* Merge back channels */
    cv::merge(chan, ycrcbImg);

    ycrcbImg.copyTo(output);
}

void inverseTransformImage(const cv::Mat& transfImg, cv::Mat& output, AxDCT_algorithm *alg){

    if(&output == nullptr) std::cerr << "\n[WARNING] Output image in inverseTransformImage function is not null. It will be erased.";

    /* Split transfImg into 3 channels */
    std::vector<cv::Mat> chan(3);
    cv::split(transfImg, chan);
    chan[0].convertTo(chan[0], CV_16S);
    chan[1].convertTo(chan[1], CV_16S);
    chan[2].convertTo(chan[2], CV_16S);

    /* Parameters for transformation */
    const int blockSize = 8;

    /********* LUMA *********/

    /* Split channel in blocks 8x8 */
    cv::Mat **tiles = splitInTiles(chan[0], blockSize);

    for(int i=0;i<chan[0].rows/blockSize;i++){
        for(int j=0;j<chan[0].cols/blockSize;j++){

            /* Dequantization step */
            y_dequantizate(tiles[i][j], tiles[i][j], alg);

            /* Do the exact IDCT */
            tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
            cv::idct(tiles[i][j], tiles[i][j]);

            /* Convert back to uint8 */
            tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[0].rows, chan[0].cols)).copyTo(chan[0]);
    chan[0].convertTo(chan[0], CV_8U);

    /**********************/

    /********* Cr *********/

    /* Split channel in blocks 8x8 */
    tiles = splitInTiles(chan[1], blockSize);

    for(int i=0;i<chan[1].rows/blockSize;i++){
        for(int j=0;j<chan[1].cols/blockSize;j++){

            /* Dequantization step */
            cr_dequantizate(tiles[i][j], tiles[i][j], alg);

            /* Do the exact IDCT */
            tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
            cv::idct(tiles[i][j], tiles[i][j]);

            /* Convert back to uint8 */
            tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[1].rows, chan[1].cols)).copyTo(chan[1]);
    chan[1].convertTo(chan[1], CV_8U);

    /**********************/

    /********* Cb *********/

    /* Split channel in blocks 8x8 */
    tiles = splitInTiles(chan[2], blockSize);

    for(int i=0;i<chan[2].rows/blockSize;i++){
        for(int j=0;j<chan[2].cols/blockSize;j++){

            /* Dequantization step */
            cb_dequantizate(tiles[i][j], tiles[i][j], alg);

            /* Do the exact IDCT */
            tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
            cv::idct(tiles[i][j], tiles[i][j]);

            /* Convert back to uint8 */
            tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[2].rows, chan[2].cols)).copyTo(chan[2]);
    chan[2].convertTo(chan[2], CV_8U);

    /**********************/

    /* Init the output matrix */
    output = cv::Mat::zeros(transfImg.rows, transfImg.cols, CV_8U);
    
    /* Merge back channels */
    cv::merge(chan, output);

    /* Converto to BGR and show image */
    cv::cvtColor(output, output, cv::COLOR_YCrCb2BGR);

}
