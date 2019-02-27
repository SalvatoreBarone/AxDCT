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

#include "core/dct.h"
#include "algorithms_list.h"
#include <cnl/fixed_point.h>
using cnl::fixed_point;

#define CHECKPOINT (std::cerr<<"\n\n\n"<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;
void test_fixed(){
    auto x = fixed_point<int16_t, -8>{1/sqrt(8)};
    std::cout<<"********************************\n\n\n";

    std::cout <<"\nto rep: " <<to_rep(x);
    std::cout <<"\nnormale: " <<x;

    std::cout<<"\n\n\n********************************";
}
int main(int argc, char** argv )
{

    // test_fixed();
    // return 0;
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    // Load img
    cv::Mat bgrImg = imread( argv[1], cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");

    // Declare an empty for dst image
    cv::Mat ycrcbImg;

    /* Convert BGR image into YCrCb */
    cv::cvtColor(bgrImg, ycrcbImg, cv::COLOR_BGR2YCrCb);

    /* Split YCbCr into 3 channels */
    std::vector<cv::Mat> chan(3);
    cv::split(ycrcbImg, chan);

    chan[0].convertTo(chan[0], CV_16S);
    // chan[1].convertTo(chan[1], CV_16S);
    // chan[2].convertTo(chan[2], CV_16S);

    /* Retrieve parameters for transformation */
    int blockSize = 8;
    // cv::Mat T, D, Q, CQ;
    // BC12::retrieveParameters(T, D, Q, CQ);

    /********* LUMA *********/

    /* Split channel in blocks 8x8 */
    cv::Mat **tiles = splitInTiles(chan[0], 8);

    for(int i=0;i<chan[0].rows/blockSize;i++){
        for(int j=0;j<chan[0].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], tiles[i][j]);
            if( i==0 && j==0) PRINT_MAT(tiles[i][j], "AxDCT");

            /* Quantization step */
            y_quantizate(tiles[i][j], tiles[i][j]);
            if( i==0 && j==0) PRINT_MAT(tiles[i][j], "y_quantizate");

            /* Dequantization step */
            y_dequantizate(tiles[i][j], tiles[i][j]);
            if( i==0 && j==0) PRINT_MAT(tiles[i][j], "y_DEquantizate");

            /* Do the exact IDCT */
            tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
            cv::idct(tiles[i][j], tiles[i][j]);
            if( i==0 && j==0) PRINT_MAT(tiles[i][j], "idct");

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
    // tiles = splitInTiles(chan[1], 8);

    // for(int i=0;i<chan[1].rows/blockSize;i++){
    //     for(int j=0;j<chan[1].cols/blockSize;j++){
            
    //         /* Do the Approximate DCT */
    //         AxDCT(tiles[i][j], tiles[i][j]);

    //         /* Quantization step */
    //         cr_quantizate(tiles[i][j], tiles[i][j]);

    //         /* Dequantization step */
    //         cr_dequantizate(tiles[i][j], tiles[i][j]);

    //         /* Do the exact IDCT */
    //         tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
    //         cv::idct(tiles[i][j], tiles[i][j]);

    //         /* Convert back to uint8 */
    //         tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
    //     }
    // }

    // /* Merge blocks 8x8 into one matrix */
    // (mergeTiles(tiles, chan[1].rows, chan[1].cols)).copyTo(chan[1]);
    // chan[1].convertTo(chan[1], CV_8U);

    // /**********************/

    // /********* Cb *********/

    // /* Split channel in blocks 8x8 */
    // tiles = splitInTiles(chan[2], 8);

    // for(int i=0;i<chan[2].rows/blockSize;i++){
    //     for(int j=0;j<chan[2].cols/blockSize;j++){
            
    //         /* Do the Approximate DCT */
    //         AxDCT(tiles[i][j], tiles[i][j]);

    //         /* Quantization step */
    //         cb_quantizate(tiles[i][j], tiles[i][j]);

    //         /* Dequantization step */
    //         cb_dequantizate(tiles[i][j], tiles[i][j]);

    //         /* Do the exact IDCT */
    //         tiles[i][j].convertTo(tiles[i][j], CV_64FC1);
    //         cv::idct(tiles[i][j], tiles[i][j]);

    //         /* Convert back to uint8 */
    //         tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
    //     }
    // }

    // /* Merge blocks 8x8 into one matrix */
    // (mergeTiles(tiles, chan[2].rows, chan[2].cols)).copyTo(chan[2]);
    // chan[2].convertTo(chan[2], CV_8U);

    /**********************/
    
    /* Merge back channels */
    cv::merge(chan, ycrcbImg);

    /* Converto to BGR and show image */
    cv::cvtColor(ycrcbImg, ycrcbImg, cv::COLOR_YCrCb2BGR);
    cv::namedWindow("Modified Image", cv::WINDOW_AUTOSIZE );
    imshow("Modified Image", ycrcbImg);
    
    cv::waitKey(0);
    return 0;
}








