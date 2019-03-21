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
 * @date   11 mar 2019
 * @brief  Implementation of mutant_eval executable function
 ******************************************************************************/

#include "core/dct.h"
#include "../utils/generic_utils.h"
#include "metrics/metrics.h"
#include "metrics/psnr_metric_eval.h"
#include "metrics/mssim_metric_eval.h"
#include "metrics/dssim_metric_eval.h"
#include "metrics/mse_metric_eval.h"
#include "metrics/ad_metric_eval.h"
#include "metrics/md_metric_eval.h"
#include "algorithms_list.h"
#include <getopt.h>
#include "nablist.h"
#include "celltype_list.h"
#include <dirent.h>
#include <climits>
#include <cstdio>
#include <unistd.h>

#define ASSIGNVAL(var, nabstr, val) \
    if(nabstr == #var) { var = val; return; }

static void usage();
static void assignGlobalValue(std::string globalvararg);
static void mapGlobalValue(const std::string&, const int);

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
    std::string globalvararg = "";
    std::string metric = "";
    int count = -1;

    double (*PEA14_compute_metric)(const cv::Mat& orig);
    double (*PEA12_compute_metric)(const cv::Mat& orig);
    double (*CB11_compute_metric)(const cv::Mat& orig);
    double (*BC12_compute_metric)(const cv::Mat& orig);
    double (*BAS08_compute_metric)(const cv::Mat& orig);
    double (*BAS09_compute_metric)(const cv::Mat& orig);
    double (*BAS11_compute_metric)(const cv::Mat& orig, const double a_param);

	while ((c = getopt(argc, argv, ":p:x:i:hlasn:f:m:r:")) != -1)
	{
		switch (c)
		{
		case 'l':
			utils::printSupportedAlgsAndMetrics();
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
            globalvararg = optarg;
            assignGlobalValue(globalvararg);
			break;

		case 'h':
			usage();
			return EXIT_SUCCESS;

        case 'f':
            folder_path = optarg;
            break;

        case 'r':
            count = atoi(optarg);
            break;

        
        case 'm':
            metric = optarg;
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

    if(( metric == "psnr")||(metric == "PSNR")){
        PEA14_compute_metric = &(metrics::PEA14_PSNR);
        PEA12_compute_metric = &(metrics::PEA12_PSNR);
        CB11_compute_metric = &(metrics::CB11_PSNR);
        BC12_compute_metric = &(metrics::BC12_PSNR);
        BAS08_compute_metric = &(metrics::BAS08_PSNR);
        BAS09_compute_metric = &(metrics::BAS09_PSNR);
        BAS11_compute_metric = &(metrics::BAS11_PSNR);
    } else if (( metric == "mse")||(metric == "MSE")){
        PEA14_compute_metric = &(metrics::PEA14_MSE);
        PEA12_compute_metric = &(metrics::PEA12_MSE);
        CB11_compute_metric = &(metrics::CB11_MSE);
        BC12_compute_metric = &(metrics::BC12_MSE);
        BAS08_compute_metric = &(metrics::BAS08_MSE);
        BAS09_compute_metric = &(metrics::BAS09_MSE);
        BAS11_compute_metric = &(metrics::BAS11_MSE);
    } else if (( metric == "ad")||(metric == "AD")){
        PEA14_compute_metric = &(metrics::PEA14_AD);
        PEA12_compute_metric = &(metrics::PEA12_AD);
        CB11_compute_metric = &(metrics::CB11_AD);
        BC12_compute_metric = &(metrics::BC12_AD);
        BAS08_compute_metric = &(metrics::BAS08_AD);
        BAS09_compute_metric = &(metrics::BAS09_AD);
        BAS11_compute_metric = &(metrics::BAS11_AD);
    } else if (( metric == "md")||(metric == "MD")){
        PEA14_compute_metric = &(metrics::PEA14_MD);
        PEA12_compute_metric = &(metrics::PEA12_MD);
        CB11_compute_metric = &(metrics::CB11_MD);
        BC12_compute_metric = &(metrics::BC12_MD);
        BAS08_compute_metric = &(metrics::BAS08_MD);
        BAS09_compute_metric = &(metrics::BAS09_MD);
        BAS11_compute_metric = &(metrics::BAS11_MD);
    } else if (( metric == "mssim")||(metric == "MSSIM")){
        PEA14_compute_metric = &(metrics::PEA14_MSSIM);
        PEA12_compute_metric = &(metrics::PEA12_MSSIM);
        CB11_compute_metric = &(metrics::CB11_MSSIM);
        BC12_compute_metric = &(metrics::BC12_MSSIM);
        BAS08_compute_metric = &(metrics::BAS08_MSSIM);
        BAS09_compute_metric = &(metrics::BAS09_MSSIM);
        BAS11_compute_metric = &(metrics::BAS11_MSSIM);
    } else if (( metric == "dssim")||(metric == "DSSIM")){
        PEA14_compute_metric = &(metrics::PEA14_DSSIM);
        PEA12_compute_metric = &(metrics::PEA12_DSSIM);
        CB11_compute_metric = &(metrics::CB11_DSSIM);
        BC12_compute_metric = &(metrics::BC12_DSSIM);
        BAS08_compute_metric = &(metrics::BAS08_DSSIM);
        BAS09_compute_metric = &(metrics::BAS09_DSSIM);
        BAS11_compute_metric = &(metrics::BAS11_DSSIM);
    } else {
        std::cout << "\nA valid metric is mandatory.";
        usage();
        return EXIT_FAILURE;
    }
    

    if(img_path != ""){

        std::vector<double> vals;

        // Load img
        cv::Mat bgrImg = imread( img_path.c_str(), cv::IMREAD_COLOR );
        assert( bgrImg.data && "No image data");

        if( algorithm == "__all") {
            vals.push_back(BC12_compute_metric(bgrImg) );
            vals.push_back(CB11_compute_metric(bgrImg) );
            vals.push_back(BAS08_compute_metric(bgrImg) );
            vals.push_back(BAS09_compute_metric(bgrImg) );
            vals.push_back(BAS11_compute_metric(bgrImg, 0.0) );
            vals.push_back(BAS11_compute_metric(bgrImg, 0.5) );
            vals.push_back(BAS11_compute_metric(bgrImg, 1.0) );
            vals.push_back(BAS11_compute_metric(bgrImg, 2.0) );
            vals.push_back(PEA12_compute_metric(bgrImg) );
            vals.push_back(PEA14_compute_metric(bgrImg) );

            utils::print_results(metric, vals, isOutputSilent);

        } else if( algorithm == "BC12" || algorithm == "bc12"){
            vals.push_back(BC12_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "CB11" || algorithm == "cb11"){
            vals.push_back(CB11_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "BAS08" || algorithm == "bas08"){
            vals.push_back(BAS08_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "BAS09" || algorithm == "bas09"){
            vals.push_back(BAS09_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "BAS11" || algorithm == "bas11"){
            vals.push_back(BAS11_compute_metric(bgrImg, a_param) );
            std::string str = std::to_string(a_param);
            str.erase (str.find_last_not_of('0') + 1, std::string::npos);
            if(a_param != 0.5) str.append("0");
            algorithm.append(" - a=" + str);
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "PEA12" || algorithm == "pea12"){
            vals.push_back(PEA12_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else if( algorithm == "PEA14" || algorithm == "pea14"){
            vals.push_back(PEA14_compute_metric(bgrImg) );
            utils::print_single_result(metric, vals, algorithm, isOutputSilent);

        } else {
            std::cout << "\nChosen algorithm (" + algorithm + ") is not supported yet.\n";
            usage();
            return EXIT_FAILURE;
        }
    } else {
        std::vector<double> vals;
        std::vector<std::string> list = utils::listFolder(folder_path);
        for (const auto& entry : list){
            std::string anImg = folder_path;
            anImg.append("/");
            anImg.append(entry);
            if( (anImg.find(".tiff") == std::string::npos) && (anImg.find(".bmp") == std::string::npos) && (anImg.find(".png") == std::string::npos) && (anImg.find(".jpg") == std::string::npos) ) continue;
            // Load img
            cv::Mat bgrImg = imread( anImg.c_str(), cv::IMREAD_COLOR );
            if( bgrImg.data == nullptr){ 
                continue;
            }

            if( algorithm == "BC12" || algorithm == "bc12"){
                vals.push_back(BC12_compute_metric(bgrImg) );

            } else if( algorithm == "CB11" || algorithm == "cb11"){
                vals.push_back(CB11_compute_metric(bgrImg) );

            } else if( algorithm == "BAS08" || algorithm == "bas08"){
                vals.push_back(BAS08_compute_metric(bgrImg) );

            } else if( algorithm == "BAS09" || algorithm == "bas09"){
                vals.push_back(BAS09_compute_metric(bgrImg) );

            } else if( algorithm == "BAS11" || algorithm == "bas11"){
                vals.push_back(BAS11_compute_metric(bgrImg, a_param) );
                std::string str = std::to_string(a_param);
                str.erase (str.find_last_not_of('0') + 1, std::string::npos);
                if(a_param != 0.5) str.append("0");
                algorithm.append(" - a=" + str);

            } else if( algorithm == "PEA12" || algorithm == "pea12"){
                vals.push_back(PEA12_compute_metric(bgrImg) );

            } else if( algorithm == "PEA14" || algorithm == "pea14"){
                vals.push_back(PEA14_compute_metric(bgrImg) );

            } else {
                std::cout << "\nChosen algorithm (" + algorithm + ") is not supported yet.\n";
                usage();
                return EXIT_FAILURE;
            }
            
            if(count-- == 0) break;

        }

        std::vector<double> mean = {0.0};
        for( auto val : vals){
            mean.at(0) += val;
        }
        mean.at(0) /= vals.size();
        utils::print_single_result(metric, mean, algorithm, isOutputSilent);

    }

    return EXIT_SUCCESS;
}

 
void usage(){
    
	std::cout << "\n\n psnr             [OPTION] [VALUE]                                \n";
	std::cout << " -i	<VALUE>		    Source image path                               \n";
	std::cout << " -f	<VALUE>		    Folder path (to run among several images)       \n";
    std::cout << " -r	<VALUE>		    Maximum number of images to evaluate within the folder       \n";
    std::cout << " -x	<VALUE>		    Chosen AxDCT algorithm                          \n";
	std::cout << " -a	    		    Compute PSNR for every algorithm                \n";
    std::cout << " -p   <VALUE>         Addition parameter for BAS11 algorithm          \n";
    std::cout << " -n   <ID>    <VALUE>	Assign a VALUE to the global variable ID        \n";
	std::cout << " -m	<VALUE>		    Metric to evaluate                              \n";
    std::cout << " -s	        	    Silent execution                        	    \n";
    std::cout << " -l	        	    List of supported algorithms and metrics        \n";
    std::cout << " -h	        	    Help                        	                \n";

	std::cout << std::endl;
}



void assignGlobalValue(std::string globalvararg){
    std::string globalIdVal [2] = { "", "" };
    std::string delimiter = " ";

    size_t pos = 0;
    int c=0;
    while ((pos = globalvararg.find(delimiter)) != std::string::npos) {
        globalIdVal[c++] = globalvararg.substr(0, pos);
        globalvararg.erase(0, pos + delimiter.length());
    }
    globalIdVal[1] = globalvararg;
    int val = std::stoi(globalIdVal[1]);
    mapGlobalValue(globalIdVal[0], val);
}

void mapGlobalValue(const std::string& id, const int val ){
    
    ASSIGNVAL(base_0, id, val)
    ASSIGNVAL(nab_0, id, val)
    ASSIGNVAL(nab_1, id, val)
    ASSIGNVAL(nab_2, id, val)
    ASSIGNVAL(nab_3, id, val)
    ASSIGNVAL(nab_4, id, val)
    ASSIGNVAL(nab_5, id, val)
    ASSIGNVAL(nab_6, id, val)
    ASSIGNVAL(nab_7, id, val)
    ASSIGNVAL(nab_8, id, val)
    ASSIGNVAL(nab_9, id, val)
    ASSIGNVAL(nab_10, id, val)
    ASSIGNVAL(nab_11, id, val)
    ASSIGNVAL(nab_12, id, val)
    ASSIGNVAL(nab_13, id, val)
    ASSIGNVAL(nab_14, id, val)
    ASSIGNVAL(nab_15, id, val)
    ASSIGNVAL(nab_16, id, val)
    ASSIGNVAL(nab_17, id, val)
    ASSIGNVAL(nab_18, id, val)
    ASSIGNVAL(nab_19, id, val)
    ASSIGNVAL(nab_20, id, val)
    ASSIGNVAL(nab_21, id, val)
    ASSIGNVAL(nab_22, id, val)
    ASSIGNVAL(nab_23, id, val)
    ASSIGNVAL(nab_24, id, val)
    ASSIGNVAL(nab_25, id, val)
    ASSIGNVAL(nab_26, id, val)
    ASSIGNVAL(nab_27, id, val)
    ASSIGNVAL(nab_28, id, val)
    ASSIGNVAL(nab_29, id, val)
    ASSIGNVAL(nab_30, id, val)
    ASSIGNVAL(nab_31, id, val)
    ASSIGNVAL(nab_32, id, val)
    ASSIGNVAL(nab_33, id, val)
    ASSIGNVAL(nab_34, id, val)
    ASSIGNVAL(nab_35, id, val)
    ASSIGNVAL(nab_36, id, val)
    ASSIGNVAL(nab_37, id, val)
    ASSIGNVAL(nab_38, id, val)
    ASSIGNVAL(nab_39, id, val)
    ASSIGNVAL(nab_40, id, val)
    ASSIGNVAL(nab_41, id, val)
    ASSIGNVAL(nab_42, id, val)
    ASSIGNVAL(nab_43, id, val)
    ASSIGNVAL(nab_44, id, val)
    ASSIGNVAL(nab_45, id, val)
    ASSIGNVAL(nab_46, id, val)
    ASSIGNVAL(nab_47, id, val)
    ASSIGNVAL(nab_48, id, val)
    ASSIGNVAL(nab_49, id, val)
    ASSIGNVAL(nab_50, id, val)
    ASSIGNVAL(nab_51, id, val)
    ASSIGNVAL(nab_52, id, val)
    ASSIGNVAL(nab_53, id, val)
    ASSIGNVAL(nab_54, id, val)
    ASSIGNVAL(nab_55, id, val)
    ASSIGNVAL(nab_56, id, val)
    ASSIGNVAL(nab_57, id, val)
    ASSIGNVAL(nab_58, id, val)
    ASSIGNVAL(nab_59, id, val)
    ASSIGNVAL(nab_60, id, val)
    ASSIGNVAL(nab_61, id, val)
    ASSIGNVAL(nab_62, id, val)
    ASSIGNVAL(nab_63, id, val)
    ASSIGNVAL(nab_64, id, val)
    ASSIGNVAL(nab_65, id, val)
    ASSIGNVAL(nab_66, id, val)
    ASSIGNVAL(nab_67, id, val)
    ASSIGNVAL(nab_68, id, val)
    ASSIGNVAL(nab_69, id, val)
    ASSIGNVAL(nab_70, id, val)
    ASSIGNVAL(nab_71, id, val)
    ASSIGNVAL(nab_72, id, val)
    ASSIGNVAL(nab_73, id, val)
    ASSIGNVAL(nab_74, id, val)
    ASSIGNVAL(nab_75, id, val)
    ASSIGNVAL(nab_76, id, val)
    ASSIGNVAL(nab_77, id, val)
    ASSIGNVAL(nab_78, id, val)
    ASSIGNVAL(nab_79, id, val)
    ASSIGNVAL(nab_80, id, val)
    ASSIGNVAL(nab_81, id, val)
    ASSIGNVAL(nab_82, id, val)
    ASSIGNVAL(nab_83, id, val)
    ASSIGNVAL(nab_84, id, val)
    ASSIGNVAL(nab_85, id, val)
    ASSIGNVAL(nab_86, id, val)
    ASSIGNVAL(nab_87, id, val)
    ASSIGNVAL(nab_88, id, val)
    ASSIGNVAL(nab_89, id, val)
    ASSIGNVAL(nab_90, id, val)
    ASSIGNVAL(nab_91, id, val)
    ASSIGNVAL(nab_92, id, val)
    ASSIGNVAL(nab_93, id, val)
    ASSIGNVAL(nab_94, id, val)
    ASSIGNVAL(nab_95, id, val)
    ASSIGNVAL(nab_96, id, val)
    ASSIGNVAL(nab_97, id, val)
    ASSIGNVAL(nab_98, id, val)
    ASSIGNVAL(nab_99, id, val)
    ASSIGNVAL(nab_100, id, val)
    ASSIGNVAL(nab_101, id, val)
    ASSIGNVAL(nab_102, id, val)
    ASSIGNVAL(nab_103, id, val)
    ASSIGNVAL(nab_104, id, val)
    ASSIGNVAL(nab_105, id, val)
    ASSIGNVAL(nab_106, id, val)
    ASSIGNVAL(nab_107, id, val)
    ASSIGNVAL(nab_108, id, val)
    ASSIGNVAL(nab_109, id, val)
    ASSIGNVAL(nab_110, id, val)
    ASSIGNVAL(nab_111, id, val)
    ASSIGNVAL(nab_112, id, val)
    ASSIGNVAL(nab_113, id, val)
    ASSIGNVAL(nab_114, id, val)
    ASSIGNVAL(nab_115, id, val)
    ASSIGNVAL(nab_116, id, val)
    ASSIGNVAL(nab_117, id, val)
    ASSIGNVAL(nab_118, id, val)
    ASSIGNVAL(nab_119, id, val)
    ASSIGNVAL(nab_120, id, val)
    ASSIGNVAL(nab_121, id, val)
    ASSIGNVAL(nab_122, id, val)
    ASSIGNVAL(nab_123, id, val)
    ASSIGNVAL(nab_124, id, val)
    ASSIGNVAL(nab_125, id, val)
    ASSIGNVAL(nab_126, id, val)
    ASSIGNVAL(nab_127, id, val)
    ASSIGNVAL(nab_128, id, val)
    ASSIGNVAL(nab_129, id, val)
    ASSIGNVAL(nab_130, id, val)
    ASSIGNVAL(nab_131, id, val)
    ASSIGNVAL(nab_132, id, val)
    ASSIGNVAL(nab_133, id, val)
    ASSIGNVAL(nab_134, id, val)
    ASSIGNVAL(nab_135, id, val)
    ASSIGNVAL(nab_136, id, val)
    ASSIGNVAL(nab_137, id, val)
    ASSIGNVAL(nab_138, id, val)
    ASSIGNVAL(cellType_0, id, val)
    ASSIGNVAL(cellType_1, id, val)
    ASSIGNVAL(cellType_2, id, val)
    ASSIGNVAL(cellType_3, id, val)
    ASSIGNVAL(cellType_4, id, val)
    ASSIGNVAL(cellType_5, id, val)
    ASSIGNVAL(cellType_6, id, val)
    ASSIGNVAL(cellType_7, id, val)
    ASSIGNVAL(cellType_8, id, val)
    ASSIGNVAL(cellType_9, id, val)
    ASSIGNVAL(cellType_10, id, val)
    ASSIGNVAL(cellType_11, id, val)
    ASSIGNVAL(cellType_12, id, val)
    ASSIGNVAL(cellType_13, id, val)
    ASSIGNVAL(cellType_14, id, val)
    ASSIGNVAL(cellType_15, id, val)
    ASSIGNVAL(cellType_16, id, val)
    ASSIGNVAL(cellType_17, id, val)
    ASSIGNVAL(cellType_18, id, val)
    ASSIGNVAL(cellType_19, id, val)
    ASSIGNVAL(cellType_20, id, val)
    ASSIGNVAL(cellType_21, id, val)
    ASSIGNVAL(cellType_22, id, val)
    ASSIGNVAL(cellType_23, id, val)
    ASSIGNVAL(cellType_24, id, val)
    ASSIGNVAL(cellType_25, id, val)
    ASSIGNVAL(cellType_26, id, val)
    ASSIGNVAL(cellType_27, id, val)
    ASSIGNVAL(cellType_28, id, val)
    ASSIGNVAL(cellType_29, id, val)
    ASSIGNVAL(cellType_30, id, val)
    ASSIGNVAL(cellType_31, id, val)
    ASSIGNVAL(cellType_32, id, val)
    ASSIGNVAL(cellType_33, id, val)
    ASSIGNVAL(cellType_34, id, val)
    ASSIGNVAL(cellType_35, id, val)
    ASSIGNVAL(cellType_36, id, val)
    ASSIGNVAL(cellType_37, id, val)
    ASSIGNVAL(cellType_38, id, val)
    ASSIGNVAL(cellType_39, id, val)
    ASSIGNVAL(cellType_40, id, val)
    ASSIGNVAL(cellType_41, id, val)
    ASSIGNVAL(cellType_42, id, val)
    ASSIGNVAL(cellType_43, id, val)
    ASSIGNVAL(cellType_44, id, val)
    ASSIGNVAL(cellType_45, id, val)
    ASSIGNVAL(cellType_46, id, val)
    ASSIGNVAL(cellType_47, id, val)
    ASSIGNVAL(cellType_48, id, val)
    ASSIGNVAL(cellType_49, id, val)
    ASSIGNVAL(cellType_50, id, val)
    ASSIGNVAL(cellType_51, id, val)
    ASSIGNVAL(cellType_52, id, val)
    ASSIGNVAL(cellType_53, id, val)
    ASSIGNVAL(cellType_54, id, val)
    ASSIGNVAL(cellType_55, id, val)
    ASSIGNVAL(cellType_56, id, val)
    ASSIGNVAL(cellType_57, id, val)
    ASSIGNVAL(cellType_58, id, val)
    ASSIGNVAL(cellType_59, id, val)
    ASSIGNVAL(cellType_60, id, val)
    ASSIGNVAL(cellType_61, id, val)
    ASSIGNVAL(cellType_62, id, val)
    ASSIGNVAL(cellType_63, id, val)
    ASSIGNVAL(cellType_64, id, val)
    ASSIGNVAL(cellType_65, id, val)
    ASSIGNVAL(cellType_66, id, val)
    ASSIGNVAL(cellType_67, id, val)
    ASSIGNVAL(cellType_68, id, val)
    ASSIGNVAL(cellType_69, id, val)
    ASSIGNVAL(cellType_70, id, val)
    ASSIGNVAL(cellType_71, id, val)
    ASSIGNVAL(cellType_72, id, val)
    ASSIGNVAL(cellType_73, id, val)
    ASSIGNVAL(cellType_74, id, val)
    ASSIGNVAL(cellType_75, id, val)
    ASSIGNVAL(cellType_76, id, val)
    ASSIGNVAL(cellType_77, id, val)
    ASSIGNVAL(cellType_78, id, val)
    ASSIGNVAL(cellType_79, id, val)
    ASSIGNVAL(cellType_80, id, val)
    ASSIGNVAL(cellType_81, id, val)
    ASSIGNVAL(cellType_82, id, val)
    ASSIGNVAL(cellType_83, id, val)
    ASSIGNVAL(cellType_84, id, val)
    ASSIGNVAL(cellType_85, id, val)
    ASSIGNVAL(cellType_86, id, val)
    ASSIGNVAL(cellType_87, id, val)
    ASSIGNVAL(cellType_88, id, val)
    ASSIGNVAL(cellType_89, id, val)
    ASSIGNVAL(cellType_90, id, val)
    ASSIGNVAL(cellType_91, id, val)
    ASSIGNVAL(cellType_92, id, val)
    ASSIGNVAL(cellType_93, id, val)
    ASSIGNVAL(cellType_94, id, val)
    ASSIGNVAL(cellType_95, id, val)
    ASSIGNVAL(cellType_96, id, val)
    ASSIGNVAL(cellType_97, id, val)
    ASSIGNVAL(cellType_98, id, val)
    ASSIGNVAL(cellType_99, id, val)
    ASSIGNVAL(cellType_100, id, val)
    ASSIGNVAL(cellType_101, id, val)
    ASSIGNVAL(cellType_102, id, val)
    ASSIGNVAL(cellType_103, id, val)
    ASSIGNVAL(cellType_104, id, val)
    ASSIGNVAL(cellType_105, id, val)
    ASSIGNVAL(cellType_106, id, val)
    ASSIGNVAL(cellType_107, id, val)
    ASSIGNVAL(cellType_108, id, val)
    ASSIGNVAL(cellType_109, id, val)
    ASSIGNVAL(cellType_110, id, val)
    ASSIGNVAL(cellType_111, id, val)
    ASSIGNVAL(cellType_112, id, val)
    ASSIGNVAL(cellType_113, id, val)
    ASSIGNVAL(cellType_114, id, val)
    ASSIGNVAL(cellType_115, id, val)
    ASSIGNVAL(cellType_116, id, val)
    ASSIGNVAL(cellType_117, id, val)
    ASSIGNVAL(cellType_118, id, val)
    ASSIGNVAL(cellType_119, id, val)
    ASSIGNVAL(cellType_120, id, val)
    ASSIGNVAL(cellType_121, id, val)
    ASSIGNVAL(cellType_122, id, val)
    ASSIGNVAL(cellType_123, id, val)
    ASSIGNVAL(cellType_124, id, val)
    ASSIGNVAL(cellType_125, id, val)
    ASSIGNVAL(cellType_126, id, val)
    ASSIGNVAL(cellType_127, id, val)
    ASSIGNVAL(cellType_128, id, val)
    ASSIGNVAL(cellType_129, id, val)
    ASSIGNVAL(cellType_130, id, val)
    ASSIGNVAL(cellType_131, id, val)
    ASSIGNVAL(cellType_132, id, val)
    ASSIGNVAL(cellType_133, id, val)
    ASSIGNVAL(cellType_134, id, val)
    ASSIGNVAL(cellType_135, id, val)
    ASSIGNVAL(cellType_136, id, val)
    ASSIGNVAL(cellType_137, id, val)
    ASSIGNVAL(cellType_138, id, val)

    std::cerr << "\nUnexpected nab value";
    assert(false);
    
}

