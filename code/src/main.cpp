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
#include <getopt.h>

#define CHECKPOINT (std::cerr<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;

void usage();
AxDCT_algorithm *stringToAlgorithm(std::string);

int main(int argc, char** argv )
{
    if( argc < 2){
        usage();
        return EXIT_FAILURE;
    }

    int c = 0;
    std::string algorithm = "";
    std::string img_path = "";

	while ((c = getopt(argc, argv, "x:i:h")) != -1)
	{
		switch (c)
		{
        case 'x':
			algorithm = optarg;
			break;

        case 'i':
			img_path = optarg;
			break;

		case 'h':
			usage();
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

    // Load img
    cv::Mat bgrImg = imread( img_path, cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");

    // Declare an empty image for transformation
    cv::Mat transfImg = bgrImg;

    // Direct and inverse transform
    AxDCT_algorithm *alg = stringToAlgorithm(algorithm);

    transformImage(bgrImg,transfImg, alg );
    inverseTransformImage(transfImg, bgrImg, alg);

    delete alg;

    // Show the approximate image 
    cv::namedWindow("Approximate Image", cv::WINDOW_AUTOSIZE );
    imshow("Approximate Image", bgrImg);
    
    cv::waitKey(0);
    return 0;
}

void usage(){
	std::cout << "\n\n axdct         [OPTION] [VALUE]                   \n";
	std::cout << " -i	<VALUE>		Source image path                   \n";
    std::cout << " -x	<VALUE>		Chosen AxDCT algorithm              \n";
    std::cout << " -h	        	Help                        	    \n";
	std::cout << std::endl;
}

AxDCT_algorithm *stringToAlgorithm(std::string algorithm){
    if( algorithm == "BC12" || algorithm == "bc12"){
        return new BC12;

    } else if( algorithm == "CB11" || algorithm == "cb11"){
        return new CB11;

    } else if( algorithm == "BAS08" || algorithm == "bas08"){
        return new BAS08;

    } else if( algorithm == "BAS09" || algorithm == "bas09"){
        return new BAS09;

    } else if( algorithm == "BAS11" || algorithm == "bas11"){
        return new BAS11;

    } else if( algorithm == "PEA12" || algorithm == "pea12"){
        return new PEA12;

    } else if( algorithm == "PEA14" || algorithm == "pea14"){
        return new PEA14;

    }
}




