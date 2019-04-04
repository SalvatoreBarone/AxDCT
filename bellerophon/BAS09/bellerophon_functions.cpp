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
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x BAS09 -m dssim";
    
    cmd.append(" -n \"nab_20 " + std::to_string(nab_20) + "\"" );
    cmd.append(" -n \"nab_21 " + std::to_string(nab_21) + "\"" );
    cmd.append(" -n \"nab_22 " + std::to_string(nab_22) + "\"" );
    cmd.append(" -n \"nab_23 " + std::to_string(nab_23) + "\"" );
    cmd.append(" -n \"nab_24 " + std::to_string(nab_24) + "\"" );
    cmd.append(" -n \"nab_25 " + std::to_string(nab_25) + "\"" );
    cmd.append(" -n \"nab_26 " + std::to_string(nab_26) + "\"" );
    cmd.append(" -n \"nab_27 " + std::to_string(nab_27) + "\"" );
    cmd.append(" -n \"nab_28 " + std::to_string(nab_28) + "\"" );
    cmd.append(" -n \"nab_29 " + std::to_string(nab_29) + "\"" );
    cmd.append(" -n \"nab_30 " + std::to_string(nab_30) + "\"" );
    cmd.append(" -n \"nab_31 " + std::to_string(nab_31) + "\"" );
    cmd.append(" -n \"nab_32 " + std::to_string(nab_32) + "\"" );
    cmd.append(" -n \"nab_33 " + std::to_string(nab_33) + "\"" );
    cmd.append(" -n \"nab_34 " + std::to_string(nab_34) + "\"" );
    cmd.append(" -n \"nab_35 " + std::to_string(nab_35) + "\"" );
    cmd.append(" -n \"nab_36 " + std::to_string(nab_36) + "\"" );
    cmd.append(" -n \"nab_37 " + std::to_string(nab_37) + "\"" );
    cmd.append(" -n \"nab_38 " + std::to_string(nab_38) + "\"" );
    cmd.append(" -n \"nab_39 " + std::to_string(nab_39) + "\"" );

    cmd.append(" -n \"cellType_20 " + std::to_string(cellType_20) + "\"" );
    cmd.append(" -n \"cellType_21 " + std::to_string(cellType_21) + "\"" );
    cmd.append(" -n \"cellType_22 " + std::to_string(cellType_22) + "\"" );
    cmd.append(" -n \"cellType_23 " + std::to_string(cellType_23) + "\"" );
    cmd.append(" -n \"cellType_24 " + std::to_string(cellType_24) + "\"" );
    cmd.append(" -n \"cellType_25 " + std::to_string(cellType_25) + "\"" );
    cmd.append(" -n \"cellType_26 " + std::to_string(cellType_26) + "\"" );
    cmd.append(" -n \"cellType_27 " + std::to_string(cellType_27) + "\"" );
    cmd.append(" -n \"cellType_28 " + std::to_string(cellType_28) + "\"" );
    cmd.append(" -n \"cellType_29 " + std::to_string(cellType_29) + "\"" );
    cmd.append(" -n \"cellType_30 " + std::to_string(cellType_30) + "\"" );
    cmd.append(" -n \"cellType_31 " + std::to_string(cellType_31) + "\"" );
    cmd.append(" -n \"cellType_32 " + std::to_string(cellType_32) + "\"" );
    cmd.append(" -n \"cellType_33 " + std::to_string(cellType_33) + "\"" );
    cmd.append(" -n \"cellType_34 " + std::to_string(cellType_34) + "\"" );
    cmd.append(" -n \"cellType_35 " + std::to_string(cellType_35) + "\"" );
    cmd.append(" -n \"cellType_36 " + std::to_string(cellType_36) + "\"" );
    cmd.append(" -n \"cellType_37 " + std::to_string(cellType_37) + "\"" );
    cmd.append(" -n \"cellType_38 " + std::to_string(cellType_38) + "\"" );
    cmd.append(" -n \"cellType_39 " + std::to_string(cellType_39) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );
    

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_20 % 10) == 0) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 1) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 2) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 3) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 4) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 5) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 6) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 7) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 8) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_20 % 10) == 9) ret += (( nab_20 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_21 % 10) == 0) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 1) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 2) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 3) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 4) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 5) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 6) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 7) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 8) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_21 % 10) == 9) ret += (( nab_21 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_22 % 10) == 0) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 1) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 2) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 3) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 4) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 5) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 6) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 7) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 8) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_22 % 10) == 9) ret += (( nab_22 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_23 % 10) == 0) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 1) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 2) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 3) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 4) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 5) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 6) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 7) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 8) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_23 % 10) == 9) ret += (( nab_23 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_24 % 10) == 0) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 1) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 2) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 3) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 4) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 5) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 6) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 7) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 8) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_24 % 10) == 9) ret += (( nab_24 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_25 % 10) == 0) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 1) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 2) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 3) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 4) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 5) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 6) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 7) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 8) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_25 % 10) == 9) ret += (( nab_25 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_26 % 10) == 0) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 1) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 2) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 3) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 4) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 5) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 6) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 7) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 8) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_26 % 10) == 9) ret += (( nab_26 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_27 % 10) == 0) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 1) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 2) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 3) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 4) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 5) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 6) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 7) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 8) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_27 % 10) == 9) ret += (( nab_27 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_28 % 10) == 0) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 1) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 2) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 3) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 4) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 5) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 6) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 7) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 8) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_28 % 10) == 9) ret += (( nab_28 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_29 % 10) == 0) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 1) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 2) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 3) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 4) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 5) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 6) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 7) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 8) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_29 % 10) == 9) ret += (( nab_29 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_30 % 10) == 0) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 1) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 2) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 3) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 4) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 5) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 6) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 7) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 8) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_30 % 10) == 9) ret += (( nab_30 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_31 % 10) == 0) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 1) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 2) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 3) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 4) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 5) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 6) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 7) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 8) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_31 % 10) == 9) ret += (( nab_31 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_32 % 10) == 0) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 1) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 2) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 3) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 4) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 5) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 6) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 7) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 8) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_32 % 10) == 9) ret += (( nab_32 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_33 % 10) == 0) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 1) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 2) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 3) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 4) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 5) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 6) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 7) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 8) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_33 % 10) == 9) ret += (( nab_33 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_34 % 10) == 0) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 1) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 2) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 3) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 4) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 5) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 6) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 7) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 8) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_34 % 10) == 9) ret += (( nab_34 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_35 % 10) == 0) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 1) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 2) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 3) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 4) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 5) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 6) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 7) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 8) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_35 % 10) == 9) ret += (( nab_35 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_36 % 10) == 0) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 1) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 2) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 3) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 4) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 5) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 6) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 7) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 8) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_36 % 10) == 9) ret += (( nab_36 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_37 % 10) == 0) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 1) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 2) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 3) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 4) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 5) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 6) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 7) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 8) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_37 % 10) == 9) ret += (( nab_37 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_38 % 10) == 0) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 1) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 2) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 3) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 4) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 5) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 6) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 7) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 8) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_38 % 10) == 9) ret += (( nab_38 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_39 % 10) == 0) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 1) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 2) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 3) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 4) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 5) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 6) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 7) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 8) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_39 % 10) == 9) ret += (( nab_39 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    

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
