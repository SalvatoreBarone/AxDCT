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
#include "algorithms_list.h"

int main(int argc, char** argv){

    assert(argc == 3 && "usage: psnr <image_path> <algorithm>\n");

    //TODO: check if algorithm exists
    //TODO: check if img exists

    std::vector<double> vals();



    return EXIT_SUCCESS;
}