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


#define CHECKPOINT (std::cerr<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

int main(int argc, char** argv )
{
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

    /* Retrieve parameters for transformation */
    int blockSize = 8;
    cv::Mat T, D, Q, CQ;
    BC12::retrieveParameters(T, D, Q, CQ);

    /********* LUMA *********/

    /* Split channel in blocks 8x8 */
    cv::Mat **tiles = splitInTiles(chan[0], 8);

    for(int i=0;i<chan[0].rows/blockSize;i++){
        for(int j=0;j<chan[0].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], T, tiles[i][j]);

            /* Quantization step */
            quantizate(tiles[i][j], D, Q, tiles[i][j]);

            /* Dequantization step */
            dequantizate(tiles[i][j], Q, tiles[i][j]);

            /* Do the exact IDCT */
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
    tiles = splitInTiles(chan[1], 8);

    for(int i=0;i<chan[1].rows/blockSize;i++){
        for(int j=0;j<chan[1].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], T, tiles[i][j]);

            /* Quantization step */
            quantizate(tiles[i][j], D, CQ, tiles[i][j]);

            /* Dequantization step */
            dequantizate(tiles[i][j], CQ, tiles[i][j]);

            /* Do the exact IDCT */
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
    tiles = splitInTiles(chan[2], 8);

    for(int i=0;i<chan[2].rows/blockSize;i++){
        for(int j=0;j<chan[2].cols/blockSize;j++){
            
            /* Do the Approximate DCT */
            AxDCT(tiles[i][j], T, tiles[i][j]);

            /* Quantization step */
            quantizate(tiles[i][j], D, CQ, tiles[i][j]);

            /* Dequantization step */
            dequantizate(tiles[i][j], CQ, tiles[i][j]);

            /* Do the exact IDCT */
            cv::idct(tiles[i][j], tiles[i][j]);

            /* Convert back to uint8 */
            tiles[i][j].convertTo(tiles[i][j], CV_8U);
            
        }
    }

    /* Merge blocks 8x8 into one matrix */
    (mergeTiles(tiles, chan[2].rows, chan[2].cols)).copyTo(chan[2]);
    chan[2].convertTo(chan[2], CV_8U);

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

template<typename T> 
void matrix_mult(const cv::Mat &A, const cv::Mat &B, cv::Mat &RES, int type){
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
    matrix_mult<int16_t>(T, tile, output, CV_16S);
    matrix_mult<int16_t>(output, T_t, output, CV_16S);
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
