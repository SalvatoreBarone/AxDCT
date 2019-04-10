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
 * @file   md_metric_eval.cpp
 * @author Andrea Aletto
 * @date   12 mar 2019
 * @brief  Implementation of md calc function
 ******************************************************************************/

#include "zybo_md_metric_eval.h"

double metrics::BC12_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BC12_transf_img = orig;
    cv::Mat BC12_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BC12_transf_img, new BC12_zybo);
    inverseTransformImage(BC12_transf_img, BC12_itransf_img, new BC12_zybo);

    // Compute PSNR and return it
    return compute_md(orig, BC12_itransf_img);
}

double metrics::CB11_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat CB11_transf_img = orig;
    cv::Mat CB11_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,CB11_transf_img, new CB11_zybo);
    inverseTransformImage(CB11_transf_img, CB11_itransf_img, new CB11_zybo);

    // Compute PSNR and return it
    return compute_md(orig, CB11_itransf_img);
}

double metrics::BAS08_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BAS08_transf_img = orig;
    cv::Mat BAS08_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS08_transf_img, new BAS08_zybo);
    inverseTransformImage(BAS08_transf_img, BAS08_itransf_img, new BAS08_zybo);

    // Compute PSNR and return it
    return compute_md(orig, BAS08_itransf_img);
}

double metrics::BAS09_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BAS09_transf_img = orig;
    cv::Mat BAS09_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS09_transf_img, new BAS09_zybo);
    inverseTransformImage(BAS09_transf_img, BAS09_itransf_img, new BAS09_zybo);

    // Compute PSNR and return it
    return compute_md(orig, BAS09_itransf_img);
}

double metrics::BAS11_zybo_MD(const cv::Mat& orig, double a_param){

    // Declare an empty image for transformation
    cv::Mat BAS11_transf_img = orig;
    cv::Mat BAS11_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS11_transf_img, new BAS11_zybo(a_param));
    inverseTransformImage(BAS11_transf_img, BAS11_itransf_img, new BAS11_zybo(a_param));

    // Compute PSNR and return it
    return compute_md(orig, BAS11_itransf_img);
}
 
double metrics::PEA12_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat PEA12_transf_img = orig;
    cv::Mat PEA12_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,PEA12_transf_img, new PEA12_zybo);
    inverseTransformImage(PEA12_transf_img, PEA12_itransf_img, new PEA12_zybo);

    // Compute PSNR and return it
    return compute_md(orig, PEA12_itransf_img);
}

double metrics::PEA14_zybo_MD(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat PEA14_transf_img = orig;
    cv::Mat PEA14_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,PEA14_transf_img, new PEA14_zybo);
    inverseTransformImage(PEA14_transf_img, PEA14_itransf_img, new PEA14_zybo);

    // Compute PSNR and return it
    return compute_md(orig, PEA14_itransf_img);
}
