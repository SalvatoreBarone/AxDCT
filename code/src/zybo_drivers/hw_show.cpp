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
 * @file   hw_show.cpp
 * @author Andrea Aletto
 * @date   4 feb 2019
 * @brief  Implementation of hw_show executable
 ******************************************************************************/

#include "core/dct.h"
#include "algorithms_list.h"
#include <getopt.h>
#include "generic_utils.h"
#include "bc12_zybo.h"
#include "cb11_zybo.h"
#include "pea12_zybo.h"
#include "pea14_zybo.h"
#include "bas08_zybo.h"
#include "bas09_zybo.h"
#include "bas11_zybo.h"

#define CHECKPOINT (std::cerr<<"\n\n\n"<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

void usage();
AxDCT_algorithm *stringToAlgorithm(std::string, bool, double a_param = -1);
void showAxDCTImage(const cv::Mat&, const std::string&, bool, double = -1, bool save = true);

int main(int argc, char** argv )
{
    if( argc == 1){
        usage();
        return EXIT_FAILURE;
    }

    int c = 0;
    std::string algorithm = "";
    double a_param = -1;
    std::string img_path = "";
    bool save = false;

	while ((c = getopt(argc, argv, ":p:x:i:hls")) != -1){
		switch (c)
		{
        case 'x':
			algorithm = optarg;
			break;

        case 'p':
            a_param = atof(optarg);
            break;

        case 'i':
			img_path = optarg;
			break;

        case 's':
			save = true;
			break;

		case 'h':
			usage();
			return EXIT_SUCCESS;
        
        case 'l':
			utils::printSupportedAlgsAndMetrics();
			return EXIT_SUCCESS;

		default:
			std::cout << "\n\nInvalid option: \n\n";
			usage();
			return EXIT_FAILURE;
		}
	}

    if( img_path == ""){
        std::cout << "\nImage path is mandatory.";
        usage();
        return EXIT_FAILURE;
    }

    if( (algorithm == "BAS11") && (a_param == -1)){
        std::cout << "\nBAS11 requires the specification of 'a' parameter. Please use -p option, too.";
        usage();
        return EXIT_FAILURE;
    }

    // Load img
    cv::Mat bgrImg = imread( img_path, cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");


    showAxDCTImage(bgrImg,algorithm, false, a_param, save);
    showAxDCTImage(bgrImg,algorithm, true, a_param, save);
    
    
    
    cv::waitKey(0);
    return 0;
}

void showAxDCTImage(const cv::Mat& bgrImg, const std::string& algorithm, bool useCustomHw, double a_param, bool save){
    
    // Declare an empty image for transformation
    cv::Mat transfImg = bgrImg;
    cv::Mat itransfImg = bgrImg;

    // Direct and inverse transform
    AxDCT_algorithm *alg = stringToAlgorithm(algorithm, useCustomHw, a_param);

    transformImage(bgrImg,transfImg, alg );
    inverseTransformImage(transfImg, itransfImg, alg);

    delete alg;

    if(save){
        imwrite( "/root/approximate_image.bmp", itransfImg );
    }
    // Show the approximate image 
    std::string winName("Approximate Image (" + algorithm );
    if(algorithm == "BAS11") {
        std::string str = std::to_string(a_param);
        str.erase (str.find_last_not_of('0') + 1, std::string::npos);
        if(a_param != 0.5) str.append("0");
        winName.append(" - a=" + str);
    }
    std::string arch = 
    winName.append(" - " + ( useCustomHw ? (std::string("software")) : (std::string("hardware")) ) );
    winName.append(")");

    cv::namedWindow(winName.c_str(), cv::WINDOW_AUTOSIZE );
    imshow(winName.c_str(), itransfImg);
}

void usage(){
	std::cout << "\n\n axdct         [OPTION] [VALUE]                                   \n";
	std::cout << " -i	<VALUE>		Source image path                                   \n";
    std::cout << " -x	<VALUE>		Chosen AxDCT algorithm                              \n";
    std::cout << " -p   <VALUE>         Addition parameter for BAS11 algorithm          \n";
    std::cout << " -l	        	    List of supported algorithms and metrics        \n";
    std::cout << " -h	        	Help                        	                    \n";
	std::cout << std::endl;

}

AxDCT_algorithm *stringToAlgorithm(std::string algorithm, bool useCustomHw, double a_param){
    if(useCustomHw) {
        if( algorithm == "BC12" || algorithm == "bc12"){
            return new BC12_zybo;

        } else if( algorithm == "CB11" || algorithm == "cb11"){
            return new CB11_zybo;

        } else if( algorithm == "BAS08" || algorithm == "bas08"){
            return new BAS08_zybo;

        } else if( algorithm == "BAS09" || algorithm == "bas09"){
            return new BAS09_zybo;

        } else if( algorithm == "BAS11" || algorithm == "bas11"){
            return new BAS11_zybo(a_param);

        } else if( algorithm == "PEA12" || algorithm == "pea12"){
            return new PEA12_zybo;

        } else if( algorithm == "PEA14" || algorithm == "pea14"){
            return new PEA14_zybo;

        } else return nullptr;
    } else {
        if( algorithm == "BC12" || algorithm == "bc12"){
            return new BC12;

        } else if( algorithm == "CB11" || algorithm == "cb11"){
            return new CB11;

        } else if( algorithm == "BAS08" || algorithm == "bas08"){
            return new BAS08;

        } else if( algorithm == "BAS09" || algorithm == "bas09"){
            return new BAS09;

        } else if( algorithm == "BAS11" || algorithm == "bas11"){
            return new BAS11(a_param);

        } else if( algorithm == "PEA12" || algorithm == "pea12"){
            return new PEA12;

        } else if( algorithm == "PEA14" || algorithm == "pea14"){
            return new PEA14;

        } else return nullptr;
    }
    
}




