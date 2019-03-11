/******************************************************************************
 * @file   bellerophon_functions.cpp
 * @author Andrea Aletto
 * @date   13 Feb 2019
 * @brief  Bellerophon functions for the DCT Algorithm Approximated through
 *         clang-chimera
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

// #include "algorithms_list.h"
// #include "algorithms/BAS08.h"
// #include "core/dct.h"
// #include "core/metrics.h"
#include <opencv2/opencv.hpp>



int nab_68 = 0;
int nab_67 = 0;
int nab_66 = 0;
int nab_65 = 0;
int nab_64 = 0;
int nab_63 = 0;
int nab_62 = 0;
int nab_61 = 0;
int nab_60 = 0;
int nab_59 = 0;
int nab_58 = 0;
int nab_57 = 0;
int nab_56 = 0;
int nab_55 = 0;

std::string exec(const char* cmd);

extern "C" double BELLERO_getError()
{
    std::string cmd = "/home/andrea/AxDCT/code/build/bin/psnr -s -i /home/andrea/lena.bmp -x BC12 -n \"nab_55 " + std::to_string(nab_55) + "\"";
    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
    // cv::Mat x = cv::Mat::zeros(8,8, CV_8U);
    return 0.0;

}


std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


