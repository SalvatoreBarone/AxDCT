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
#include "nablist.h"
#include <dirent.h>
#include <climits>
#include <cstdio>
#include <unistd.h>

#define ASSIGNVAL(var, nabstr, val) \
    if(nabstr == #var) { var = val; return; }

static void usage();
static void printSupportedAlgs();
static void print_results(std::vector<double>, bool silent = false);
static void print_single_result(std::vector<double>, std::string, bool silent = false);
static void assignNabValue(std::string nabarg);
static void mapNabValue(const std::string& nabid, const int nabVal );
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
    std::string folder_path = "";
    bool isOutputSilent = false;
    std::string nabarg = "";

	while ((c = getopt(argc, argv, "p:x:i::hlasn:f:")) != -1)
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
    mapNabValue(nabIdVal[0], nabval);
}

void mapNabValue(const std::string& nabid, const int nabVal ){
    
    ASSIGNVAL(nab_0, nabid, nabVal)
    ASSIGNVAL(nab_1, nabid, nabVal)
    ASSIGNVAL(nab_2, nabid, nabVal)
    ASSIGNVAL(nab_3, nabid, nabVal)
    ASSIGNVAL(nab_4, nabid, nabVal)
    ASSIGNVAL(nab_5, nabid, nabVal)
    ASSIGNVAL(nab_6, nabid, nabVal)
    ASSIGNVAL(nab_7, nabid, nabVal)
    ASSIGNVAL(nab_8, nabid, nabVal)
    ASSIGNVAL(nab_9, nabid, nabVal)
    ASSIGNVAL(nab_10, nabid, nabVal)
    ASSIGNVAL(nab_11, nabid, nabVal)
    ASSIGNVAL(nab_12, nabid, nabVal)
    ASSIGNVAL(nab_13, nabid, nabVal)
    ASSIGNVAL(nab_14, nabid, nabVal)
    ASSIGNVAL(nab_15, nabid, nabVal)
    ASSIGNVAL(nab_16, nabid, nabVal)
    ASSIGNVAL(nab_17, nabid, nabVal)
    ASSIGNVAL(nab_18, nabid, nabVal)
    ASSIGNVAL(nab_19, nabid, nabVal)
    ASSIGNVAL(nab_20, nabid, nabVal)
    ASSIGNVAL(nab_21, nabid, nabVal)
    ASSIGNVAL(nab_22, nabid, nabVal)
    ASSIGNVAL(nab_23, nabid, nabVal)
    ASSIGNVAL(nab_24, nabid, nabVal)
    ASSIGNVAL(nab_25, nabid, nabVal)
    ASSIGNVAL(nab_26, nabid, nabVal)
    ASSIGNVAL(nab_27, nabid, nabVal)
    ASSIGNVAL(nab_28, nabid, nabVal)
    ASSIGNVAL(nab_29, nabid, nabVal)
    ASSIGNVAL(nab_30, nabid, nabVal)
    ASSIGNVAL(nab_31, nabid, nabVal)
    ASSIGNVAL(nab_32, nabid, nabVal)
    ASSIGNVAL(nab_33, nabid, nabVal)
    ASSIGNVAL(nab_34, nabid, nabVal)
    ASSIGNVAL(nab_35, nabid, nabVal)
    ASSIGNVAL(nab_36, nabid, nabVal)
    ASSIGNVAL(nab_37, nabid, nabVal)
    ASSIGNVAL(nab_38, nabid, nabVal)
    ASSIGNVAL(nab_39, nabid, nabVal)
    ASSIGNVAL(nab_40, nabid, nabVal)
    ASSIGNVAL(nab_41, nabid, nabVal)
    ASSIGNVAL(nab_42, nabid, nabVal)
    ASSIGNVAL(nab_43, nabid, nabVal)
    ASSIGNVAL(nab_44, nabid, nabVal)
    ASSIGNVAL(nab_45, nabid, nabVal)
    ASSIGNVAL(nab_46, nabid, nabVal)
    ASSIGNVAL(nab_47, nabid, nabVal)
    ASSIGNVAL(nab_48, nabid, nabVal)
    ASSIGNVAL(nab_49, nabid, nabVal)
    ASSIGNVAL(nab_50, nabid, nabVal)
    ASSIGNVAL(nab_51, nabid, nabVal)
    ASSIGNVAL(nab_52, nabid, nabVal)
    ASSIGNVAL(nab_53, nabid, nabVal)
    ASSIGNVAL(nab_54, nabid, nabVal)
    ASSIGNVAL(nab_55, nabid, nabVal)
    ASSIGNVAL(nab_56, nabid, nabVal)
    ASSIGNVAL(nab_57, nabid, nabVal)
    ASSIGNVAL(nab_58, nabid, nabVal)
    ASSIGNVAL(nab_59, nabid, nabVal)
    ASSIGNVAL(nab_60, nabid, nabVal)
    ASSIGNVAL(nab_61, nabid, nabVal)
    ASSIGNVAL(nab_62, nabid, nabVal)
    ASSIGNVAL(nab_63, nabid, nabVal)
    ASSIGNVAL(nab_64, nabid, nabVal)
    ASSIGNVAL(nab_65, nabid, nabVal)
    ASSIGNVAL(nab_66, nabid, nabVal)
    ASSIGNVAL(nab_67, nabid, nabVal)
    ASSIGNVAL(nab_68, nabid, nabVal)
    ASSIGNVAL(nab_69, nabid, nabVal)
    ASSIGNVAL(nab_70, nabid, nabVal)
    ASSIGNVAL(nab_71, nabid, nabVal)
    ASSIGNVAL(nab_72, nabid, nabVal)
    ASSIGNVAL(nab_73, nabid, nabVal)
    ASSIGNVAL(nab_74, nabid, nabVal)
    ASSIGNVAL(nab_75, nabid, nabVal)
    ASSIGNVAL(nab_76, nabid, nabVal)
    ASSIGNVAL(nab_77, nabid, nabVal)
    ASSIGNVAL(nab_78, nabid, nabVal)
    ASSIGNVAL(nab_79, nabid, nabVal)
    ASSIGNVAL(nab_80, nabid, nabVal)
    ASSIGNVAL(nab_81, nabid, nabVal)
    ASSIGNVAL(nab_82, nabid, nabVal)
    ASSIGNVAL(nab_83, nabid, nabVal)
    ASSIGNVAL(nab_84, nabid, nabVal)
    ASSIGNVAL(nab_85, nabid, nabVal)
    ASSIGNVAL(nab_86, nabid, nabVal)
    ASSIGNVAL(nab_87, nabid, nabVal)
    ASSIGNVAL(nab_88, nabid, nabVal)
    ASSIGNVAL(nab_89, nabid, nabVal)
    ASSIGNVAL(nab_90, nabid, nabVal)
    ASSIGNVAL(nab_91, nabid, nabVal)
    ASSIGNVAL(nab_92, nabid, nabVal)
    ASSIGNVAL(nab_93, nabid, nabVal)
    ASSIGNVAL(nab_94, nabid, nabVal)
    ASSIGNVAL(nab_95, nabid, nabVal)
    ASSIGNVAL(nab_96, nabid, nabVal)
    ASSIGNVAL(nab_97, nabid, nabVal)
    ASSIGNVAL(nab_98, nabid, nabVal)
    ASSIGNVAL(nab_99, nabid, nabVal)
    ASSIGNVAL(nab_100, nabid, nabVal)
    ASSIGNVAL(nab_101, nabid, nabVal)
    ASSIGNVAL(nab_102, nabid, nabVal)
    ASSIGNVAL(nab_103, nabid, nabVal)
    ASSIGNVAL(nab_104, nabid, nabVal)
    ASSIGNVAL(nab_105, nabid, nabVal)
    ASSIGNVAL(nab_106, nabid, nabVal)
    ASSIGNVAL(nab_107, nabid, nabVal)
    ASSIGNVAL(nab_108, nabid, nabVal)
    ASSIGNVAL(nab_109, nabid, nabVal)
    ASSIGNVAL(nab_110, nabid, nabVal)
    ASSIGNVAL(nab_111, nabid, nabVal)
    ASSIGNVAL(nab_112, nabid, nabVal)
    ASSIGNVAL(nab_113, nabid, nabVal)
    ASSIGNVAL(nab_114, nabid, nabVal)
    ASSIGNVAL(nab_115, nabid, nabVal)
    ASSIGNVAL(nab_116, nabid, nabVal)
    ASSIGNVAL(nab_117, nabid, nabVal)
    ASSIGNVAL(nab_118, nabid, nabVal)
    ASSIGNVAL(nab_119, nabid, nabVal)
    ASSIGNVAL(nab_120, nabid, nabVal)
    ASSIGNVAL(nab_121, nabid, nabVal)
    ASSIGNVAL(nab_122, nabid, nabVal)
    ASSIGNVAL(nab_123, nabid, nabVal)
    ASSIGNVAL(nab_124, nabid, nabVal)
    ASSIGNVAL(nab_125, nabid, nabVal)
    ASSIGNVAL(nab_126, nabid, nabVal)
    ASSIGNVAL(nab_127, nabid, nabVal)
    ASSIGNVAL(nab_128, nabid, nabVal)

    std::cerr << "\nUnexpected nab value";
    assert(false);
    
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