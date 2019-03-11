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
#include <dirent.h>
#include <climits>
#include <cstdio>
#include <unistd.h>

void usage();
void printSupportedAlgs();
void print_results(std::vector<double>, bool silent = false);
void print_single_result(std::vector<double>, std::string, bool silent = false);
static std::vector<std::string> listFolder(std::string dir);

double BC12_PSNR(const cv::Mat& orig);
double CB11_PSNR(const cv::Mat& orig);
double BAS08_PSNR(const cv::Mat& orig);
double BAS09_PSNR(const cv::Mat& orig);
double BAS11_PSNR(const cv::Mat& orig, double a_param);
double PEA12_PSNR(const cv::Mat& orig);
double PEA14_PSNR(const cv::Mat& orig);

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
    std::string folder_path = "";

	while ((c = getopt(argc, argv, "p:x:i:hlasf:")) != -1)
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

		case 'h':
			usage();
			return EXIT_SUCCESS;

        case 'f':
            folder_path = optarg;
            break;

		default:
			std::cout << "\n\nInvalid option: \n\n";
			usage();
			return EXIT_FAILURE;
		}
	}

    if((img_path == "") && (folder_path == "")){
        std::cout << "\nImage path or folder path is mandatory.";
        usage();
        return EXIT_FAILURE;
    }

    if( algorithm == "BAS11" && a_param == -1){
        std::cout << "\nBAS11 requires the specification of 'a' parameter. Please use -p option, too.";
        usage();
        return EXIT_FAILURE;
    }

    if(img_path != ""){

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
    } else {
        std::vector<double> vals;
        std::vector<std::string> list = listFolder(folder_path);
        for (const auto& entry : list){
            std::string anImg = folder_path;
            anImg.append("/");
            anImg.append(entry);
            if(anImg.find(".bmp") == std::string::npos) continue;
            // Load img
            cv::Mat bgrImg = imread( anImg.c_str(), cv::IMREAD_COLOR );
            if( bgrImg.data == nullptr){ 
                continue;
            }

            if( algorithm == "BC12" || algorithm == "bc12"){
                vals.push_back(BC12_PSNR(bgrImg) );

            } else if( algorithm == "CB11" || algorithm == "cb11"){
                vals.push_back(CB11_PSNR(bgrImg) );

            } else if( algorithm == "BAS08" || algorithm == "bas08"){
                vals.push_back(BAS08_PSNR(bgrImg) );

            } else if( algorithm == "BAS09" || algorithm == "bas09"){
                vals.push_back(BAS09_PSNR(bgrImg) );

            } else if( algorithm == "BAS11" || algorithm == "bas11"){
                vals.push_back(BAS11_PSNR(bgrImg, a_param) );
                std::string str = std::to_string(a_param);
                str.erase (str.find_last_not_of('0') + 1, std::string::npos);
                if(a_param != 0.5) str.append("0");
                algorithm.append(" - a=" + str);

            } else if( algorithm == "PEA12" || algorithm == "pea12"){
                vals.push_back(PEA12_PSNR(bgrImg) );

            } else if( algorithm == "PEA14" || algorithm == "pea14"){
                vals.push_back(PEA14_PSNR(bgrImg) );

            } else {
                std::cout << "\nChosen algorithm (" + algorithm + ") is not supported yet.\n";
                usage();
                return EXIT_FAILURE;
            }

        }
        std::vector<double> mean = {0.0}; //media dei psnr
        for( auto val : vals){
            mean.at(0) += val;
        }
        mean.at(0) /= vals.size();
        print_single_result(mean, algorithm, isOutputSilent);

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

std::vector<std::string> listFolder(std::string path){
    std::vector<std::string>ret;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (path.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string elem = ent->d_name;
            ret.push_back(elem);
        }
        closedir (dir);
        
    }
    return ret;
}