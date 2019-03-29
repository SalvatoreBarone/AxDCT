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
 * Internal functions declaration
*/


static double getPSNR(const cv::Mat& I1, const cv::Mat& I2);
static double getMSE(const cv::Mat& I1, const cv::Mat& I2);
static double compute_ad(const cv::Mat& orig, const cv::Mat& target, int component);
static double compute_md(const cv::Mat& orig, const cv::Mat& target, int component);
static cv::Scalar getMSSIM( const cv::Mat& i1, const cv::Mat& i2);

/**
 * Public functions implementation
*/

double compute_mse(const cv::Mat& orig, const cv::Mat& target){
    return getMSE(orig, target);
}

double compute_psnr(const cv::Mat& orig, const cv::Mat& target){
    return getPSNR(orig,target);
}

double compute_ad(const cv::Mat& orig, const cv::Mat& target){
    double ad[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) ad[i] = compute_ad(orig, target, i);
    
    double ret = ( ad[0] + ad[1] + ad[2] )/3;
    
    return ret;
    return 0.0;

}

double compute_md(const cv::Mat& orig, const cv::Mat& target){
    double md[3] = {0.0, 0.0, 0.0};

    for( int i=0; i<3; i++ ) md[i] = compute_md(orig, target, i);
    
    if(md[0] > md[1]){
        if(md[2] > md[0] ) return md[2];
        else return md[0];
    } else {
        if(md[2] > md[1] ) return md[2];
        else return md[1];
    }

}

double compute_reduction(const double exact_param, const double inexact_param, const int nab, const int n_bit){
    double v_ext = n_bit * exact_param;
    double v_inxt = nab*inexact_param + (n_bit - nab)*exact_param;

    return (v_inxt-v_ext)/v_ext * 100;
}

double compute_mssim(const cv::Mat& orig, const cv::Mat& target){
    cv::Scalar mssim = getMSSIM(orig,target);
    // std::cout << "\nmssim components:\n" <<mssim.val[0] <<"\n" << mssim.val[1] << "\n" << mssim.val[2] << "\n\n\n";
    return ((mssim.val[0] + mssim.val[1] + mssim.val[2])/3);
}

double compute_dssim(const cv::Mat& orig, const cv::Mat& target){
    return (1-compute_mssim(orig,target));
}


/**
 * Internal functions implementation
*/

static double getPSNR(const cv::Mat& I1, const cv::Mat& I2){
    return ( 10.0*log10((255*255)/getMSE(I1, I2)) );
}

static double getMSE(const cv::Mat& I1, const cv::Mat& I2){
    cv::Mat s1;
    cv::absdiff(I1, I2, s1);       // |I1 - I2|
    s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
    s1 = s1.mul(s1);           // |I1 - I2|^2

    cv::Scalar s = sum(s1);         // sum elements per channel

    double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels

    if( sse <= 1e-10)
        return 0;
    else
        return ( sse /(double)(I1.channels() * I1.total()) );
    
}

static double compute_ad(const cv::Mat& orig, const cv::Mat& target, int component){
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

static double compute_md(const cv::Mat& orig, const cv::Mat& target, int component){
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

static cv::Scalar getMSSIM( const cv::Mat& i1, const cv::Mat& i2){
    const double C1 = 6.5025, C2 = 58.5225;
    /***************************** INITS **********************************/
    int d     = CV_32F;

    cv::Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);

    cv::Mat I2_2   = I2.mul(I2);        // I2^2
    cv::Mat I1_2   = I1.mul(I1);        // I1^2
    cv::Mat I1_I2  = I1.mul(I2);        // I1 * I2

    /*************************** END INITS **********************************/

    cv::Mat mu1, mu2;   // PRELIMINARY COMPUTING
    cv::GaussianBlur(I1, mu1, cv::Size(11, 11), 1.5);
    cv::GaussianBlur(I2, mu2, cv::Size(11, 11), 1.5);

    cv::Mat mu1_2   =   mu1.mul(mu1);
    cv::Mat mu2_2   =   mu2.mul(mu2);
    cv::Mat mu1_mu2 =   mu1.mul(mu2);

    cv::Mat sigma1_2, sigma2_2, sigma12;

    GaussianBlur(I1_2, sigma1_2, cv::Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;

    GaussianBlur(I2_2, sigma2_2, cv::Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;

    GaussianBlur(I1_I2, sigma12, cv::Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;

    ///////////////////////////////// FORMULA ////////////////////////////////
    cv::Mat t1, t2, t3;

    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

    cv::Mat ssim_map;
    divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;

    cv::Scalar mssim = mean( ssim_map ); // mssim = average of ssim map
    return mssim;
}
