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

#ifndef _GENERIC_UTILS_H
#define _GENERIC_UTILS_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include <dirent.h>
#include <climits>
#include <cstdio>
#include <unistd.h>

namespace utils {

    extern std::vector<std::string> supported_algorithms;
    extern std::vector<std::string> supported_metrics;

    void printSupportedAlgsAndMetrics();
    void print_results(std::string metric, std::vector<double>, bool silent = false);
    void print_single_result(std::string metric, std::vector<double>, std::string, bool silent = false);
    std::vector<std::string> listFolder(std::string dir);

} /* end of namespace utils */

#endif /* _GENERIC_UTILS_H */