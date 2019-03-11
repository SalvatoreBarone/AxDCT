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
 * @file   psnr.cpp
 * @author Andrea Aletto
 * @date   28 feb 2019
 * @brief  Implementation of psnr executable function
 ******************************************************************************/

#include "core/dct.h"
#include "core/metrics.h"
#include "algorithms_list.h"
#include <getopt.h>

void usage();
void printSupportedAlgs();
void print_results(std::vector<double>, bool silent = false);
void print_single_result(std::vector<double>, std::string, bool silent = false);
void assignNabValue(std::string nabarg);

double BC12_PSNR(const cv::Mat& orig);
double CB11_PSNR(const cv::Mat& orig);
double BAS08_PSNR(const cv::Mat& orig);
double BAS09_PSNR(const cv::Mat& orig);
double BAS11_PSNR(const cv::Mat& orig, double a_param);
double PEA12_PSNR(const cv::Mat& orig);
double PEA14_PSNR(const cv::Mat& orig);


extern int nab_68;
extern int nab_67;
extern int nab_66;
extern int nab_65;
extern int nab_64;
extern int nab_63;
extern int nab_62;
extern int nab_61;
extern int nab_60;
extern int nab_59;
extern int nab_58;
extern int nab_57;
extern int nab_56;
extern int nab_55;

static std::vector<std::string> supported_algorithms = {"BC12\t\t", "CB11\t\t", "BAS08\t", "BAS09\t", "BAS11 (a=0.0)", "BAS11 (a=0.5)", "BAS11 (a=1.0)", "BAS11 (a=2.0)", "PEA12\t", "PEA14\t" };

int main(int argc, char** argv){
    if( argc < 2){
        usage();
        return EXIT_FAILURE;
    }

    int c = 0;
    std::string algorithm = "";
    double a_param = -1;
    std::string img_path = "";
    bool isOutputSilent = false;
    std::string nabarg = "";

	while ((c = getopt(argc, argv, "p:x:i:hlasn:")) != -1)
	{
		switch (c)
		{
		case 'l':
			printSupportedAlgs();
			return EXIT_SUCCESS;

		case 'a':
			algorithm = "__all";
			break;
        case 's':
            isOutputSilent = true;
            break;

        case 'x':
			algorithm = optarg;
			break;

        case 'p':
            a_param = atof(optarg);
            break;

        case 'i':
			img_path = optarg;
			break;
        case 'n':
            nabarg = optarg;
            assignNabValue(nabarg);
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

    if( algorithm == "BAS11" && a_param == -1){
        std::cout << "\nBAS11 requires the specification of 'a' parameter. Please use -p option, too.";
        usage();
        return EXIT_FAILURE;
    }

    std::vector<double> vals;

    // Load img
    cv::Mat bgrImg = imread( img_path.c_str(), cv::IMREAD_COLOR );
    assert( bgrImg.data && "No image data");

    if( algorithm == "__all") {
        vals.push_back(BC12_PSNR(bgrImg) );
        vals.push_back(CB11_PSNR(bgrImg) );
        vals.push_back(BAS08_PSNR(bgrImg) );
        vals.push_back(BAS09_PSNR(bgrImg) );
        vals.push_back(BAS11_PSNR(bgrImg, 0.0) );
        vals.push_back(BAS11_PSNR(bgrImg, 0.5) );
        vals.push_back(BAS11_PSNR(bgrImg, 1.0) );
        vals.push_back(BAS11_PSNR(bgrImg, 2.0) );
        vals.push_back(PEA12_PSNR(bgrImg) );
        vals.push_back(PEA14_PSNR(bgrImg) );

        print_results(vals, isOutputSilent);

    } else if( algorithm == "BC12" || algorithm == "bc12"){
        vals.push_back(BC12_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "CB11" || algorithm == "cb11"){
        vals.push_back(CB11_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "BAS08" || algorithm == "bas08"){
        vals.push_back(BAS08_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "BAS09" || algorithm == "bas09"){
        vals.push_back(BAS09_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "BAS11" || algorithm == "bas11"){
        vals.push_back(BAS11_PSNR(bgrImg, a_param) );
        std::string str = std::to_string(a_param);
        str.erase (str.find_last_not_of('0') + 1, std::string::npos);
        if(a_param != 0.5) str.append("0");
        algorithm.append(" - a=" + str);
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "PEA12" || algorithm == "pea12"){
        vals.push_back(PEA12_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else if( algorithm == "PEA14" || algorithm == "pea14"){
        vals.push_back(PEA14_PSNR(bgrImg) );
        print_single_result(vals, algorithm, isOutputSilent);

    } else {
        std::cout << "\nChosen algorithm (" + algorithm + ") is not supported yet.\n";
        usage();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

double BC12_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BC12_transf_img = orig;
    cv::Mat BC12_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BC12_transf_img, new BC12 );
    inverseTransformImage(BC12_transf_img, BC12_itransf_img, new BC12);

    // Compute PSNR and return it
    return compute_psnr(orig, BC12_itransf_img);
}

double CB11_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat CB11_transf_img = orig;
    cv::Mat CB11_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,CB11_transf_img, new CB11 );
    inverseTransformImage(CB11_transf_img, CB11_itransf_img, new CB11);

    // Compute PSNR and return it
    return compute_psnr(orig, CB11_itransf_img);
}

double BAS08_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BAS08_transf_img = orig;
    cv::Mat BAS08_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS08_transf_img, new BAS08 );
    inverseTransformImage(BAS08_transf_img, BAS08_itransf_img, new BAS08);

    // Compute PSNR and return it
    return compute_psnr(orig, BAS08_itransf_img);
}

double BAS09_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat BAS09_transf_img = orig;
    cv::Mat BAS09_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS09_transf_img, new BAS09 );
    inverseTransformImage(BAS09_transf_img, BAS09_itransf_img, new BAS09);

    // Compute PSNR and return it
    return compute_psnr(orig, BAS09_itransf_img);
}

double BAS11_PSNR(const cv::Mat& orig, double a_param){

    // Declare an empty image for transformation
    cv::Mat BAS11_transf_img = orig;
    cv::Mat BAS11_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,BAS11_transf_img, new BAS11(a_param) );
    inverseTransformImage(BAS11_transf_img, BAS11_itransf_img, new BAS11(a_param) );

    // Compute PSNR and return it
    return compute_psnr(orig, BAS11_itransf_img);
}
 
double PEA12_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat PEA12_transf_img = orig;
    cv::Mat PEA12_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,PEA12_transf_img, new PEA12 );
    inverseTransformImage(PEA12_transf_img, PEA12_itransf_img, new PEA12);

    // Compute PSNR and return it
    return compute_psnr(orig, PEA12_itransf_img);
}

double PEA14_PSNR(const cv::Mat& orig){

    // Declare an empty image for transformation
    cv::Mat PEA14_transf_img = orig;
    cv::Mat PEA14_itransf_img = orig;

    // Direct and inverse transform
    transformImage(orig,PEA14_transf_img, new PEA14 );
    inverseTransformImage(PEA14_transf_img, PEA14_itransf_img, new PEA14);

    // Compute PSNR and return it
    return compute_psnr(orig, PEA14_itransf_img);
}

void print_results(std::vector<double> vals, bool silent){
    if(silent) return;
    std::cout << "\n\n************** PSNR **************\n\n";

    for(int i=0; i<vals.size(); i++){
        std::cout << "   " << supported_algorithms.at(i) << "\t" << vals.at(i) <<std::endl;
    }

    std::cout << "\n**********************************\n\n";
}

void print_single_result(std::vector<double> vals, std::string algorithm, bool silent){
    if(silent) {
        std::cout << vals.at(0);
    } else {

        for (std::string::size_type i=0; i<algorithm.length(); i++) algorithm[i]=std::toupper(algorithm[i], *(new std::locale) );

        std::cout << "\n\n************** PSNR **************\n\n";
        std::cout << "   " << algorithm << "\t" << vals.at(0) <<std::endl;
        std::cout << "\n**********************************\n\n";
    }
}
 
void usage(){
	std::cout << "\n\n psnr         [OPTION] [VALUE]                    \n";
	std::cout << " -i	<VALUE>		Source image path                   \n";
    std::cout << " -x	<VALUE>		Chosen AxDCT algorithm              \n";
	std::cout << " -a	    		Compute PSNR for every algorithm    \n";
    std::cout << " -l	        	List of supported algorithms	    \n";
    std::cout << " -h	        	Help                        	    \n";
	std::cout << std::endl;
}

void printSupportedAlgs(){
    std::cout << "\nSupported algorithms:\n";
    for (auto alg : supported_algorithms){
        std::cout <<"\n- " << alg; 
    }
}

void assignNabValue(std::string nabarg){
    std::string nabIdVal [2] = { "", "" };
    std::string delimiter = " ";

    size_t pos = 0;
    int c=0;
    while ((pos = nabarg.find(delimiter)) != std::string::npos) {
        nabIdVal[c++] = nabarg.substr(0, pos);
        nabarg.erase(0, pos + delimiter.length());
    }
    nabIdVal[1] = nabarg;
    int nabval = std::stoi(nabIdVal[1]);
    if(nabIdVal[0] == "nab_55") {
        nab_55 = nabval;
    } else {
        std::cerr << "\nUnexpected nab value";
        assert(false);
    }
}
