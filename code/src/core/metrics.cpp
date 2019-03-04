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

/**
 * METRICS for one component
*/

double compute_mse(const cv::Mat& orig, const cv::Mat& target, int component){
    assert(((component == 0) || (component == 1) || (component == 2)) && "component must be 0, 1 or 2");

    double ret = 0.0;

    cv::Mat origCh[3];
    cv::Mat targetCh[3];

    cv::split(orig, origCh);
    cv::split(target, targetCh);

    if( component == 0){
        origCh[1].deallocate();
        origCh[2].deallocate();
        targetCh[1].deallocate();
        targetCh[2].deallocate();
    } else if (component == 1){
        origCh[0].deallocate();
        origCh[2].deallocate();
        targetCh[0].deallocate();
        targetCh[2].deallocate();
    } else {
        origCh[1].deallocate();
        origCh[0].deallocate();
        targetCh[1].deallocate();
        targetCh[0].deallocate();
    }
    

    cv::Mat origDouble;
    origCh[component].convertTo(origDouble, CV_64FC1);   
    
    cv::Mat targetDouble;   
    targetCh[component].convertTo(targetDouble, CV_64FC1);

    for(int row=0; row<orig.rows; row++){
        for(int col=0; col<orig.cols; col++){
            ret += ( (origDouble.at<double>(row,col) - targetDouble.at<double>(row,col))*(origDouble.at<double>(row,col) - targetDouble.at<double>(row,col)) );
        }
    }

    origDouble.deallocate();
    targetDouble.deallocate();
    origCh[component].deallocate();
    targetCh[component].deallocate();
    
    ret /= (orig.rows * orig.cols);
    return ret;
}

double compute_psnr(const cv::Mat& orig, const cv::Mat& target, int component){
    double mse = compute_mse(orig, target, component);
    double psnr = 10*log10( 255*255/mse );
    // std::cout << "\n\npsnr for component" <<component <<" is: --> " <<psnr <<" <--\n"; 
    return psnr;
}

double compute_ad(const cv::Mat& orig, const cv::Mat& target, int component){
    assert(((component == 0) || (component == 1) || (component == 2)) && "component must be 0, 1 or 2");

    double ret = 0.0;

    cv::Mat origCh[3];
    cv::Mat targetCh[3];

    cv::split(orig, origCh);
    cv::split(target, targetCh);

    if( component == 0){
        origCh[1].deallocate();
        origCh[2].deallocate();
        targetCh[1].deallocate();
        targetCh[2].deallocate();
    } else if (component == 1){
        origCh[0].deallocate();
        origCh[2].deallocate();
        targetCh[0].deallocate();
        targetCh[2].deallocate();
    } else {
        origCh[1].deallocate();
        origCh[0].deallocate();
        targetCh[1].deallocate();
        targetCh[0].deallocate();
    }
    

    cv::Mat origDouble;
    origCh[component].convertTo(origDouble, CV_64FC1);   
    
    cv::Mat targetDouble;   
    targetCh[component].convertTo(targetDouble, CV_64FC1);

    for(int row=0; row<orig.rows; row++){
        for(int col=0; col<orig.cols; col++){
            ret += ( origDouble.at<double>(row,col) - targetDouble.at<double>(row,col) );
        }
    }

    origDouble.deallocate();
    targetDouble.deallocate();
    origCh[component].deallocate();
    targetCh[component].deallocate();
    
    ret /= (orig.rows * orig.cols);
    return ret;
}

double compute_md(const cv::Mat& orig, const cv::Mat& target, int component){
    assert(((component == 0) || (component == 1) || (component == 2)) && "component must be 0, 1 or 2");

    double ret = 0.0;

    cv::Mat origCh[3];
    cv::Mat targetCh[3];

    cv::split(orig, origCh);
    cv::split(target, targetCh);

    if( component == 0){
        origCh[1].deallocate();
        origCh[2].deallocate();
        targetCh[1].deallocate();
        targetCh[2].deallocate();
    } else if (component == 1){
        origCh[0].deallocate();
        origCh[2].deallocate();
        targetCh[0].deallocate();
        targetCh[2].deallocate();
    } else {
        origCh[1].deallocate();
        origCh[0].deallocate();
        targetCh[1].deallocate();
        targetCh[0].deallocate();
    }

    for(int row=0; row<orig.rows; row++){
        for(int col=0; col<orig.cols; col++){
            double diff = (double)(origCh[component].at<unsigned char>(row,col) - targetCh[component].at<unsigned char>(row,col));
            if( ret < diff) ret = diff;
        }
    }

    origCh[component].deallocate();
    targetCh[component].deallocate();
    
    return ret;
}

/**
 * METRICS for the whole image
*/

double compute_mse(const cv::Mat& orig, const cv::Mat& target){
    double mse[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) mse[i] = compute_mse(orig, target, i);
    
    double ret = ( mse[0] + mse[1] + mse[2] )/3;
    
    return ret;
}

double compute_psnr(const cv::Mat& orig, const cv::Mat& target){
    
    double psnr[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) psnr[i] = compute_psnr(orig, target, i);
    
    double ret = ( psnr[0] + psnr[1] + psnr[2] )/3;
    
    return ret;
}

double compute_ad(const cv::Mat& orig, const cv::Mat& target){
    double ad[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) ad[i] = compute_ad(orig, target, i);
    
    double ret = ( ad[0] + ad[1] + ad[2] )/3;
    
    return ret;
}

double compute_md(const cv::Mat& orig, const cv::Mat& target){
    double md[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) md[i] = compute_md(orig, target, i);
    
    double max = -1;

    if(md[0] > md[1]){
        if(md[2] > md[0] ) return md[2];
        else return md[0];
    } else {
        if(md[2] > md[1] ) return md[2];
        else return md[1];
    }

}

/**
 * Other METRICS
*/

double compute_reduction(const double exact_param, const double inexact_param, const int nab, const int n_bit){
    double v_ext = n_bit * exact_param;
    double v_inxt = nab*inexact_param + (n_bit - nab)*exact_param;

    return (v_inxt-v_ext)/v_ext * 100;
}


