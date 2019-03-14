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
 * @file   generic_utils.h
 * @author Andrea Aletto
 * @date   11 mar 2019
 * @brief  Declaration of generic utils functions
 ******************************************************************************/

#include "generic_utils.h"

namespace utils {
    std::vector<std::string> supported_algorithms = {"BC12\t\t", "CB11\t\t", "BAS08\t", "BAS09\t", "BAS11 (a=0.0)", "BAS11 (a=0.5)", "BAS11 (a=1.0)", "BAS11 (a=2.0)", "PEA12\t", "PEA14\t" };    
    std::vector<std::string> supported_metrics = {"PSNR", "MSE", "MD", "AD", "MSSIM" };
}

void utils::printSupportedAlgsAndMetrics(){
    std::cout << "\nSupported algorithms:\n";
    for (auto alg : supported_algorithms){
        std::cout <<"\n- " << alg; 
    }

    std::cout << "\n\nSupported metrics:\n";
    for (auto m : supported_metrics){
        std::cout <<"\n- " << m; 
    }
}

std::vector<std::string> utils::listFolder(std::string path){
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

void utils::print_results(std::string metric, std::vector<double> vals, bool silent){
    if(silent) return;
    
    std::cout << "\n\n************** " <<std::toupper(metric, *(new std::locale)) <<" **************\n\n";

    for(int i=0; i<vals.size(); i++){
        std::cout << "   " << supported_algorithms.at(i) << "\t" << vals.at(i) <<std::endl;
    }

    std::cout << "\n**********************************\n\n";
}

void utils::print_single_result(std::string metric, std::vector<double> vals, std::string algorithm, bool silent){
    if(silent) {
        std::cout << vals.at(0);
    } else {

        for (std::string::size_type i=0; i<algorithm.length(); i++) algorithm[i]=std::toupper(algorithm[i], *(new std::locale) );

        std::cout << "\n\n************** " <<metric <<" **************\n\n";
        std::cout << "   " << algorithm << "\t" << vals.at(0) <<std::endl;
        std::cout << "\n**********************************\n\n";
    }
}