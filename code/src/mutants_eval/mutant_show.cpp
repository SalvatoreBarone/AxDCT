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
 * @file   mutant_show.cpp
 * @author Andrea Aletto
 * @date   4 feb 2019
 * @brief  Implementation of mutant_show executable
 ******************************************************************************/

#include "../core/dct.h"
#include "../algorithms_list.h"
#include <getopt.h>
#include "generic_utils.h"
#include "nablist.h"
#include "celltype_list.h"

#define CHECKPOINT (std::cerr<<"\n\n\n"<<__PRETTY_FUNCTION__<<__LINE__<<std::endl);
#define PRINT_MAT(mat, msg) std::cout<< std::endl <<msg <<":" <<std::endl <<mat <<std::endl;
#define ASSIGNVAL(var, nabstr, val) \
    if(nabstr == #var) { var = val; return; }

void usage();
AxDCT_algorithm *stringToAlgorithm(std::string, double a_param = -1);
void showAxDCTImage(const cv::Mat&, const std::string&, double = -1);

static void assignGlobalValue(std::string globalvararg);
static void mapGlobalValue(const std::string&, const int);

int main(int argc, char** argv )
{
    if( argc == 1){
        usage();
        return EXIT_FAILURE;
    }

    int c = 0;
    std::string algorithm = "";
    double a_param = -1;
    std::string globalvararg = "";
    std::string img_path = "";

	while ((c = getopt(argc, argv, ":p:x:i:hlan:")) != -1){
		switch (c)
		{
        case 'x':
			algorithm = optarg;
			break;

        case 'a':
			algorithm = "__all";
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
        
        case 'l':
			utils::printSupportedAlgsAndMetrics();
			return EXIT_SUCCESS;

        case 'n':
            globalvararg = optarg;
            assignGlobalValue(globalvararg);
			break;

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

    if( algorithm == "__all") {
        showAxDCTImage(bgrImg,"BC12");
        showAxDCTImage(bgrImg,"CB11");
        showAxDCTImage(bgrImg,"BAS08");
        showAxDCTImage(bgrImg,"BAS09");
        showAxDCTImage(bgrImg,"BAS11", 0);
        showAxDCTImage(bgrImg,"BAS11", 0.5);
        showAxDCTImage(bgrImg,"BAS11", 1.0);
        showAxDCTImage(bgrImg,"BAS11", 2.0);
        showAxDCTImage(bgrImg,"PEA12");
        showAxDCTImage(bgrImg,"PEA14");
    } else {
        showAxDCTImage(bgrImg,algorithm,a_param);
    }
    
    
    cv::waitKey(0);
    return 0;
}

void showAxDCTImage(const cv::Mat& bgrImg, const std::string& algorithm, double a_param){
    
    // Declare an empty image for transformation
    cv::Mat transfImg = bgrImg;
    cv::Mat itransfImg = bgrImg;

    // Direct and inverse transform
    AxDCT_algorithm *alg = stringToAlgorithm(algorithm,a_param);

    transformImage(bgrImg,transfImg, alg );
    inverseTransformImage(transfImg, itransfImg, alg);

    delete alg;

    // Show the approximate image 
    std::string winName("Approximate Image (" + algorithm );
    if(algorithm == "BAS11") {
        std::string str = std::to_string(a_param);
        str.erase (str.find_last_not_of('0') + 1, std::string::npos);
        if(a_param != 0.5) str.append("0");
        winName.append(" - a=" + str);
    }
    winName.append(")");

    cv::namedWindow(winName.c_str(), cv::WINDOW_AUTOSIZE );
    imshow(winName.c_str(), itransfImg);
}

void usage(){
	std::cout << "\n\n axdct         [OPTION] [VALUE]                                   \n";
	std::cout << " -i	<VALUE>		Source image path                                   \n";
    std::cout << " -x	<VALUE>		Chosen AxDCT algorithm                              \n";
    std::cout << " -n   <ID>    <VALUE>	Assign a VALUE to the global variable ID        \n";
    std::cout << " -p   <VALUE>         Addition parameter for BAS11 algorithm          \n";
    std::cout << " -a	        	Compute AxDCT image for every supported algorithm   \n";
    std::cout << " -l	        	    List of supported algorithms and metrics        \n";
    std::cout << " -h	        	Help                        	                    \n";
	std::cout << std::endl;

}

AxDCT_algorithm *stringToAlgorithm(std::string algorithm, double a_param){
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
    ASSIGNVAL(cellType_0, id, val%10)
    ASSIGNVAL(cellType_1, id, val%10)
    ASSIGNVAL(cellType_2, id, val%10)
    ASSIGNVAL(cellType_3, id, val%10)
    ASSIGNVAL(cellType_4, id, val%10)
    ASSIGNVAL(cellType_5, id, val%10)
    ASSIGNVAL(cellType_6, id, val%10)
    ASSIGNVAL(cellType_7, id, val%10)
    ASSIGNVAL(cellType_8, id, val%10)
    ASSIGNVAL(cellType_9, id, val%10)
    ASSIGNVAL(cellType_10, id, val%10)
    ASSIGNVAL(cellType_11, id, val%10)
    ASSIGNVAL(cellType_12, id, val%10)
    ASSIGNVAL(cellType_13, id, val%10)
    ASSIGNVAL(cellType_14, id, val%10)
    ASSIGNVAL(cellType_15, id, val%10)
    ASSIGNVAL(cellType_16, id, val%10)
    ASSIGNVAL(cellType_17, id, val%10)
    ASSIGNVAL(cellType_18, id, val%10)
    ASSIGNVAL(cellType_19, id, val%10)
    ASSIGNVAL(cellType_20, id, val%10)
    ASSIGNVAL(cellType_21, id, val%10)
    ASSIGNVAL(cellType_22, id, val%10)
    ASSIGNVAL(cellType_23, id, val%10)
    ASSIGNVAL(cellType_24, id, val%10)
    ASSIGNVAL(cellType_25, id, val%10)
    ASSIGNVAL(cellType_26, id, val%10)
    ASSIGNVAL(cellType_27, id, val%10)
    ASSIGNVAL(cellType_28, id, val%10)
    ASSIGNVAL(cellType_29, id, val%10)
    ASSIGNVAL(cellType_30, id, val%10)
    ASSIGNVAL(cellType_31, id, val%10)
    ASSIGNVAL(cellType_32, id, val%10)
    ASSIGNVAL(cellType_33, id, val%10)
    ASSIGNVAL(cellType_34, id, val%10)
    ASSIGNVAL(cellType_35, id, val%10)
    ASSIGNVAL(cellType_36, id, val%10)
    ASSIGNVAL(cellType_37, id, val%10)
    ASSIGNVAL(cellType_38, id, val%10)
    ASSIGNVAL(cellType_39, id, val%10)
    ASSIGNVAL(cellType_40, id, val%10)
    ASSIGNVAL(cellType_41, id, val%10)
    ASSIGNVAL(cellType_42, id, val%10)
    ASSIGNVAL(cellType_43, id, val%10)
    ASSIGNVAL(cellType_44, id, val%10)
    ASSIGNVAL(cellType_45, id, val%10)
    ASSIGNVAL(cellType_46, id, val%10)
    ASSIGNVAL(cellType_47, id, val%10)
    ASSIGNVAL(cellType_48, id, val%10)
    ASSIGNVAL(cellType_49, id, val%10)
    ASSIGNVAL(cellType_50, id, val%10)
    ASSIGNVAL(cellType_51, id, val%10)
    ASSIGNVAL(cellType_52, id, val%10)
    ASSIGNVAL(cellType_53, id, val%10)
    ASSIGNVAL(cellType_54, id, val%10)
    ASSIGNVAL(cellType_55, id, val%10)
    ASSIGNVAL(cellType_56, id, val%10)
    ASSIGNVAL(cellType_57, id, val%10)
    ASSIGNVAL(cellType_58, id, val%10)
    ASSIGNVAL(cellType_59, id, val%10)
    ASSIGNVAL(cellType_60, id, val%10)
    ASSIGNVAL(cellType_61, id, val%10)
    ASSIGNVAL(cellType_62, id, val%10)
    ASSIGNVAL(cellType_63, id, val%10)
    ASSIGNVAL(cellType_64, id, val%10)
    ASSIGNVAL(cellType_65, id, val%10)
    ASSIGNVAL(cellType_66, id, val%10)
    ASSIGNVAL(cellType_67, id, val%10)
    ASSIGNVAL(cellType_68, id, val%10)
    ASSIGNVAL(cellType_69, id, val%10)
    ASSIGNVAL(cellType_70, id, val%10)
    ASSIGNVAL(cellType_71, id, val%10)
    ASSIGNVAL(cellType_72, id, val%10)
    ASSIGNVAL(cellType_73, id, val%10)
    ASSIGNVAL(cellType_74, id, val%10)
    ASSIGNVAL(cellType_75, id, val%10)
    ASSIGNVAL(cellType_76, id, val%10)
    ASSIGNVAL(cellType_77, id, val%10)
    ASSIGNVAL(cellType_78, id, val%10)
    ASSIGNVAL(cellType_79, id, val%10)
    ASSIGNVAL(cellType_80, id, val%10)
    ASSIGNVAL(cellType_81, id, val%10)
    ASSIGNVAL(cellType_82, id, val%10)
    ASSIGNVAL(cellType_83, id, val%10)
    ASSIGNVAL(cellType_84, id, val%10)
    ASSIGNVAL(cellType_85, id, val%10)
    ASSIGNVAL(cellType_86, id, val%10)
    ASSIGNVAL(cellType_87, id, val%10)
    ASSIGNVAL(cellType_88, id, val%10)
    ASSIGNVAL(cellType_89, id, val%10)
    ASSIGNVAL(cellType_90, id, val%10)
    ASSIGNVAL(cellType_91, id, val%10)
    ASSIGNVAL(cellType_92, id, val%10)
    ASSIGNVAL(cellType_93, id, val%10)
    ASSIGNVAL(cellType_94, id, val%10)
    ASSIGNVAL(cellType_95, id, val%10)
    ASSIGNVAL(cellType_96, id, val%10)
    ASSIGNVAL(cellType_97, id, val%10)
    ASSIGNVAL(cellType_98, id, val%10)
    ASSIGNVAL(cellType_99, id, val%10)
    ASSIGNVAL(cellType_100, id, val%10)
    ASSIGNVAL(cellType_101, id, val%10)
    ASSIGNVAL(cellType_102, id, val%10)
    ASSIGNVAL(cellType_103, id, val%10)
    ASSIGNVAL(cellType_104, id, val%10)
    ASSIGNVAL(cellType_105, id, val%10)
    ASSIGNVAL(cellType_106, id, val%10)
    ASSIGNVAL(cellType_107, id, val%10)
    ASSIGNVAL(cellType_108, id, val%10)
    ASSIGNVAL(cellType_109, id, val%10)
    ASSIGNVAL(cellType_110, id, val%10)
    ASSIGNVAL(cellType_111, id, val%10)
    ASSIGNVAL(cellType_112, id, val%10)
    ASSIGNVAL(cellType_113, id, val%10)
    ASSIGNVAL(cellType_114, id, val%10)
    ASSIGNVAL(cellType_115, id, val%10)
    ASSIGNVAL(cellType_116, id, val%10)
    ASSIGNVAL(cellType_117, id, val%10)
    ASSIGNVAL(cellType_118, id, val%10)
    ASSIGNVAL(cellType_119, id, val%10)
    ASSIGNVAL(cellType_120, id, val%10)
    ASSIGNVAL(cellType_121, id, val%10)
    ASSIGNVAL(cellType_122, id, val%10)
    ASSIGNVAL(cellType_123, id, val%10)
    ASSIGNVAL(cellType_124, id, val%10)
    ASSIGNVAL(cellType_125, id, val%10)
    ASSIGNVAL(cellType_126, id, val%10)
    ASSIGNVAL(cellType_127, id, val%10)
    ASSIGNVAL(cellType_128, id, val%10)
    ASSIGNVAL(cellType_129, id, val%10)
    ASSIGNVAL(cellType_130, id, val%10)
    ASSIGNVAL(cellType_131, id, val%10)
    ASSIGNVAL(cellType_132, id, val%10)
    ASSIGNVAL(cellType_133, id, val%10)
    ASSIGNVAL(cellType_134, id, val%10)
    ASSIGNVAL(cellType_135, id, val%10)
    ASSIGNVAL(cellType_136, id, val%10)
    ASSIGNVAL(cellType_137, id, val%10)
    ASSIGNVAL(cellType_138, id, val%10)

    std::cerr << "\nUnexpected nab value";
    assert(false);
    
}


