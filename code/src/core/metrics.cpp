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
 * @file   metrics.cpp
 * @author Andrea Aletto
 * @date   28 feb 2019
 * @brief  Implementation of evaluating metrics for transformed images
 ******************************************************************************/
#include "metrics.h"

double compute_psnr(const cv::Mat& orig, const cv::Mat& target){
    
    double error[3] = {0.0, 0.0, 0.0};
    
    cv::Mat origCh[3];
    cv::Mat targetCh[3];

    cv::split(orig, origCh);
    cv::split(target, targetCh);

    for(int i=0; i<3; i++){
        cv::Mat origDouble;
        origCh[i].convertTo(origDouble, CV_64FC1);   
        
        cv::Mat targetDouble;   
        targetCh[i].convertTo(targetDouble, CV_64FC1);

        for(int row=0; row<orig.rows; row++){
            for(int col=0; col<orig.cols; col++){
                error[i] += ( (origDouble.at<double>(row,col) - targetDouble.at<double>(row,col))*(origDouble.at<double>(row,col) - targetDouble.at<double>(row,col)) );
            }
        }

        origDouble.deallocate();
        targetDouble.deallocate();
        origCh[i].deallocate();
        targetCh[i].deallocate();
    }
    
    double mse = ( error[0] + error[1] + error[2] );
    mse /= (orig.cols * orig.rows);

    double psnr = 10*log10( 255*255/mse );
    
    return psnr;
}