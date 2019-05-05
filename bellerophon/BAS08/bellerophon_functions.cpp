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
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x BAS08 -m dssim";
    
    cmd.append(" -n \"nab_0 " + std::to_string(nab_0) + "\"" );
    cmd.append(" -n \"nab_1 " + std::to_string(nab_1) + "\"" );
    cmd.append(" -n \"nab_2 " + std::to_string(nab_2) + "\"" );
    cmd.append(" -n \"nab_3 " + std::to_string(nab_3) + "\"" );
    cmd.append(" -n \"nab_4 " + std::to_string(nab_4) + "\"" );
    cmd.append(" -n \"nab_5 " + std::to_string(nab_5) + "\"" );
    cmd.append(" -n \"nab_6 " + std::to_string(nab_6) + "\"" );
    cmd.append(" -n \"nab_7 " + std::to_string(nab_7) + "\"" );
    cmd.append(" -n \"nab_8 " + std::to_string(nab_8) + "\"" );
    cmd.append(" -n \"nab_9 " + std::to_string(nab_9) + "\"" );
    cmd.append(" -n \"nab_10 " + std::to_string(nab_10) + "\"" );
    cmd.append(" -n \"nab_11 " + std::to_string(nab_11) + "\"" );
    cmd.append(" -n \"nab_12 " + std::to_string(nab_12) + "\"" );
    cmd.append(" -n \"nab_13 " + std::to_string(nab_13) + "\"" );
    cmd.append(" -n \"nab_14 " + std::to_string(nab_14) + "\"" );
    cmd.append(" -n \"nab_15 " + std::to_string(nab_15) + "\"" );
    cmd.append(" -n \"nab_16 " + std::to_string(nab_16) + "\"" );
    cmd.append(" -n \"nab_17 " + std::to_string(nab_17) + "\"" );
    cmd.append(" -n \"nab_18 " + std::to_string(nab_18) + "\"" );
    cmd.append(" -n \"nab_19 " + std::to_string(nab_19) + "\"" );

    cmd.append(" -n \"cellType_0 " + std::to_string(cellType_0) + "\"" );
    cmd.append(" -n \"cellType_1 " + std::to_string(cellType_1) + "\"" );
    cmd.append(" -n \"cellType_2 " + std::to_string(cellType_2) + "\"" );
    cmd.append(" -n \"cellType_3 " + std::to_string(cellType_3) + "\"" );
    cmd.append(" -n \"cellType_4 " + std::to_string(cellType_4) + "\"" );
    cmd.append(" -n \"cellType_5 " + std::to_string(cellType_5) + "\"" );
    cmd.append(" -n \"cellType_6 " + std::to_string(cellType_6) + "\"" );
    cmd.append(" -n \"cellType_7 " + std::to_string(cellType_7) + "\"" );
    cmd.append(" -n \"cellType_8 " + std::to_string(cellType_8) + "\"" );
    cmd.append(" -n \"cellType_9 " + std::to_string(cellType_9) + "\"" );
    cmd.append(" -n \"cellType_10 " + std::to_string(cellType_10) + "\"" );
    cmd.append(" -n \"cellType_11 " + std::to_string(cellType_11) + "\"" );
    cmd.append(" -n \"cellType_12 " + std::to_string(cellType_12) + "\"" );
    cmd.append(" -n \"cellType_13 " + std::to_string(cellType_13) + "\"" );
    cmd.append(" -n \"cellType_14 " + std::to_string(cellType_14) + "\"" );
    cmd.append(" -n \"cellType_15 " + std::to_string(cellType_15) + "\"" );
    cmd.append(" -n \"cellType_16 " + std::to_string(cellType_16) + "\"" );
    cmd.append(" -n \"cellType_17 " + std::to_string(cellType_17) + "\"" );
    cmd.append(" -n \"cellType_18 " + std::to_string(cellType_18) + "\"" );
    cmd.append(" -n \"cellType_19 " + std::to_string(cellType_19) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );
    

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_0 % 10) == 0) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 1) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 2) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 3) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 4) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 5) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 6) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 7) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 8) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_0 % 10) == 9) ret += (( nab_0 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_1 % 10) == 0) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 1) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 2) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 3) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 4) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 5) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 6) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 7) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 8) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_1 % 10) == 9) ret += (( nab_1 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_2 % 10) == 0) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 1) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 2) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 3) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 4) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 5) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 6) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 7) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 8) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_2 % 10) == 9) ret += (( nab_2 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_3 % 10) == 0) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 1) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 2) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 3) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 4) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 5) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 6) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 7) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 8) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_3 % 10) == 9) ret += (( nab_3 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_4 % 10) == 0) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 1) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 2) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 3) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 4) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 5) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 6) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 7) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 8) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_4 % 10) == 9) ret += (( nab_4 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_5 % 10) == 0) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 1) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 2) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 3) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 4) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 5) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 6) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 7) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 8) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_5 % 10) == 9) ret += (( nab_5 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_6 % 10) == 0) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 1) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 2) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 3) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 4) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 5) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 6) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 7) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 8) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_6 % 10) == 9) ret += (( nab_6 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_7 % 10) == 0) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 1) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 2) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 3) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 4) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 5) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 6) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 7) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 8) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_7 % 10) == 9) ret += (( nab_7 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_8 % 10) == 0) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 1) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 2) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 3) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 4) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 5) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 6) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 7) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 8) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_8 % 10) == 9) ret += (( nab_8 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_9 % 10) == 0) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 1) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 2) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 3) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 4) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 5) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 6) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 7) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 8) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_9 % 10) == 9) ret += (( nab_9 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_10 % 10) == 0) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 1) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 2) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 3) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 4) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 5) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 6) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 7) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 8) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_10 % 10) == 9) ret += (( nab_10 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_11 % 10) == 0) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 1) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 2) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 3) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 4) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 5) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 6) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 7) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 8) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_11 % 10) == 9) ret += (( nab_11 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_12 % 10) == 0) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 1) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 2) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 3) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 4) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 5) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 6) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 7) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 8) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_12 % 10) == 9) ret += (( nab_12 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_13 % 10) == 0) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 1) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 2) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 3) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 4) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 5) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 6) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 7) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 8) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_13 % 10) == 9) ret += (( nab_13 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_14 % 10) == 0) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 1) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 2) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 3) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 4) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 5) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 6) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 7) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 8) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_14 % 10) == 9) ret += (( nab_14 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_15 % 10) == 0) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 1) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 2) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 3) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 4) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 5) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 6) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 7) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 8) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_15 % 10) == 9) ret += (( nab_15 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_16 % 10) == 0) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 1) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 2) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 3) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 4) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 5) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 6) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 7) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 8) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_16 % 10) == 9) ret += (( nab_16 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_17 % 10) == 0) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 1) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 2) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 3) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 4) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 5) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 6) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 7) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 8) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_17 % 10) == 9) ret += (( nab_17 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_18 % 10) == 0) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 1) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 2) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 3) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 4) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 5) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 6) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 7) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 8) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_18 % 10) == 9) ret += (( nab_18 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_19 % 10) == 0) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 1) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 2) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 3) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 4) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 5) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 6) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 7) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 8) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_19 % 10) == 9) ret += (( nab_19 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    

    ret /= (20.0 * 16 * FA_TRANSISTOR_COUNT);
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
