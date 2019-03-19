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
 * @file   bellerophon_functions.cpp
 * @author Andrea Aletto
 * @date   11 mar 2019
 * @brief  Implementation of Bellerophon's fitness functions
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

#include "../bellero_nablist.h"
#include "../bellero_celltype_list.h"

#include <inexact_adders.h>
using namespace inexact_adders::metrics;

static std::string exec(const char* cmd);

extern "C" double BELLERO_getError()
{
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x BC12 -m dssim";
    
    cmd.append(" -n \"nab_55 " + std::to_string(nab_55) + "\"" );
    cmd.append(" -n \"nab_56 " + std::to_string(nab_56) + "\"" );
    cmd.append(" -n \"nab_57 " + std::to_string(nab_57) + "\"" );
    cmd.append(" -n \"nab_58 " + std::to_string(nab_58) + "\"" );
    cmd.append(" -n \"nab_59 " + std::to_string(nab_59) + "\"" );
    cmd.append(" -n \"nab_60 " + std::to_string(nab_60) + "\"" );
    cmd.append(" -n \"nab_61 " + std::to_string(nab_61) + "\"" );
    cmd.append(" -n \"nab_62 " + std::to_string(nab_62) + "\"" );
    cmd.append(" -n \"nab_63 " + std::to_string(nab_63) + "\"" );
    cmd.append(" -n \"nab_64 " + std::to_string(nab_64) + "\"" );
    cmd.append(" -n \"nab_65 " + std::to_string(nab_65) + "\"" );
    cmd.append(" -n \"nab_66 " + std::to_string(nab_66) + "\"" );
    cmd.append(" -n \"nab_67 " + std::to_string(nab_67) + "\"" );
    cmd.append(" -n \"nab_68 " + std::to_string(nab_68) + "\"" );
    cmd.append(" -n \"cellType_55 " + std::to_string(cellType_55) + "\"" );
    cmd.append(" -n \"cellType_56 " + std::to_string(cellType_56) + "\"" );
    cmd.append(" -n \"cellType_57 " + std::to_string(cellType_57) + "\"" );
    cmd.append(" -n \"cellType_58 " + std::to_string(cellType_58) + "\"" );
    cmd.append(" -n \"cellType_59 " + std::to_string(cellType_59) + "\"" );
    cmd.append(" -n \"cellType_60 " + std::to_string(cellType_60) + "\"" );
    cmd.append(" -n \"cellType_61 " + std::to_string(cellType_61) + "\"" );
    cmd.append(" -n \"cellType_62 " + std::to_string(cellType_62) + "\"" );
    cmd.append(" -n \"cellType_63 " + std::to_string(cellType_63) + "\"" );
    cmd.append(" -n \"cellType_64 " + std::to_string(cellType_64) + "\"" );
    cmd.append(" -n \"cellType_65 " + std::to_string(cellType_65) + "\"" );
    cmd.append(" -n \"cellType_66 " + std::to_string(cellType_66) + "\"" );
    cmd.append(" -n \"cellType_67 " + std::to_string(cellType_67) + "\"" );
    cmd.append(" -n \"cellType_68 " + std::to_string(cellType_68) + "\"" );
    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );
    

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_55 % 10) == 0) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 1) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 2) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 3) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 4) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 5) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 6) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 7) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 8) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_55 % 10) == 9) ret += (( nab_55 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_56 % 10) == 0) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 1) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 2) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 3) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 4) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 5) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 6) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 7) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 8) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_56 % 10) == 9) ret += (( nab_56 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_57 % 10) == 0) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 1) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 2) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 3) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 4) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 5) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 6) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 7) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 8) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_57 % 10) == 9) ret += (( nab_57 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_58 % 10) == 0) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 1) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 2) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 3) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 4) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 5) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 6) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 7) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 8) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_58 % 10) == 9) ret += (( nab_58 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_59 % 10) == 0) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 1) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 2) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 3) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 4) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 5) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 6) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 7) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 8) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_59 % 10) == 9) ret += (( nab_59 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_60 % 10) == 0) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 1) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 2) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 3) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 4) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 5) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 6) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 7) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 8) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_60 % 10) == 9) ret += (( nab_60 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_61 % 10) == 0) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 1) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 2) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 3) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 4) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 5) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 6) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 7) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 8) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_61 % 10) == 9) ret += (( nab_61 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_62 % 10) == 0) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 1) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 2) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 3) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 4) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 5) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 6) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 7) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 8) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_62 % 10) == 9) ret += (( nab_62 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_63 % 10) == 0) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 1) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 2) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 3) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 4) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 5) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 6) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 7) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 8) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_63 % 10) == 9) ret += (( nab_63 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_64 % 10) == 0) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 1) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 2) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 3) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 4) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 5) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 6) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 7) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 8) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_64 % 10) == 9) ret += (( nab_64 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_65 % 10) == 0) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 1) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 2) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 3) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 4) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 5) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 6) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 7) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 8) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_65 % 10) == 9) ret += (( nab_65 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_66 % 10) == 0) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 1) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 2) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 3) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 4) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 5) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 6) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 7) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 8) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_66 % 10) == 9) ret += (( nab_66 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_67 % 10) == 0) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 1) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 2) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 3) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 4) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 5) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 6) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 7) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 8) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_67 % 10) == 9) ret += (( nab_67 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_68 % 10) == 0) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 1) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 2) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 3) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 4) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 5) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 6) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 7) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 8) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_68 % 10) == 9) ret += (( nab_68 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    ret /= (14.0 * 16 * FA_TRANSISTOR_COUNT);
    return ret;
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


