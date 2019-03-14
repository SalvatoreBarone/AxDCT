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
 * @file   image_eval.cpp
 * @author Andrea Aletto
 * @date   28 feb 2019
 * @brief  Implementation of image_eval executable
 ******************************************************************************/

#include "core/dct.h"
#include "metrics/metrics.h"
#include "metrics/psnr_metric_eval.h"
#include "metrics/mse_metric_eval.h"
#include "metrics/ad_metric_eval.h"
#include "metrics/md_metric_eval.h"
#include "metrics/mssim_metric_eval.h"
#include "utils/generic_utils.h"
#include "algorithms_list.h"
#include <getopt.h>

void usage();


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
    std::string metric = "";

    double (*PEA14_compute_metric)(const cv::Mat& orig);
    double (*PEA12_compute_metric)(const cv::Mat& orig);
    double (*CB11_compute_metric)(const cv::Mat& orig);
    double (*BC12_compute_metric)(const cv::Mat& orig);
    double (*BAS08_compute_metric)(const cv::Mat& orig);
    double (*BAS09_compute_metric)(const cv::Mat& orig);
    double (*BAS11_compute_metric)(const cv::Mat& orig, const double a_param);

	while ((c = getopt(argc, argv, ":p:x:i:hlasf:m:")) != -1)
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

		case 'h':
			usage();
			return EXIT_SUCCESS;

        case 'f':
            folder_path = optarg;
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
            if(anImg.find(".bmp") == std::string::npos) continue;
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

        }
        std::vector<double> mean = {0.0}; //media dei psnr
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
    std::cout << " -x	<VALUE>		    Chosen AxDCT algorithm                          \n";
	std::cout << " -a	    		    Compute PSNR for every algorithm                \n";
    std::cout << " -p   <VALUE>         Addition parameter for BAS11 algorithm          \n";
	std::cout << " -m	<VALUE>		    Metric to evaluate                              \n";
    std::cout << " -s	        	    Silent execution                        	    \n";
    std::cout << " -l	        	    List of supported algorithms and metrics        \n";
    std::cout << " -h	        	    Help                        	                \n";

	std::cout << std::endl;
}

