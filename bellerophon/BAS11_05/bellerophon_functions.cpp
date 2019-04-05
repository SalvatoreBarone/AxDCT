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
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x BAS11 -a 0.5 -m dssim";
    
    cmd.append(" -n \"nab_58 " + std::to_string(nab_58) + "\"" );
    cmd.append(" -n \"nab_57 " + std::to_string(nab_57) + "\"" );
    cmd.append(" -n \"nab_56 " + std::to_string(nab_56) + "\"" );
    cmd.append(" -n \"nab_55 " + std::to_string(nab_55) + "\"" );
    cmd.append(" -n \"nab_54 " + std::to_string(nab_54) + "\"" );
    cmd.append(" -n \"nab_53 " + std::to_string(nab_53) + "\"" );
    cmd.append(" -n \"nab_52 " + std::to_string(nab_52) + "\"" );
    cmd.append(" -n \"nab_51 " + std::to_string(nab_51) + "\"" );
    cmd.append(" -n \"nab_50 " + std::to_string(nab_50) + "\"" );
    cmd.append(" -n \"nab_49 " + std::to_string(nab_49) + "\"" );
    cmd.append(" -n \"nab_48 " + std::to_string(nab_48) + "\"" );
    cmd.append(" -n \"nab_47 " + std::to_string(nab_47) + "\"" );
    cmd.append(" -n \"nab_46 " + std::to_string(nab_46) + "\"" );
    cmd.append(" -n \"nab_45 " + std::to_string(nab_45) + "\"" );
    cmd.append(" -n \"nab_44 " + std::to_string(nab_44) + "\"" );
    cmd.append(" -n \"nab_43 " + std::to_string(nab_43) + "\"" );
    cmd.append(" -n \"nab_42 " + std::to_string(nab_42) + "\"" );
    cmd.append(" -n \"nab_41 " + std::to_string(nab_41) + "\"" );
    cmd.append(" -n \"nab_40 " + std::to_string(nab_40) + "\"" );

    cmd.append(" -n \"cellType_58 " + std::to_string(cellType_58) + "\"" );
    cmd.append(" -n \"cellType_57 " + std::to_string(cellType_57) + "\"" );
    cmd.append(" -n \"cellType_56 " + std::to_string(cellType_56) + "\"" );
    cmd.append(" -n \"cellType_55 " + std::to_string(cellType_55) + "\"" );
    cmd.append(" -n \"cellType_54 " + std::to_string(cellType_54) + "\"" );
    cmd.append(" -n \"cellType_53 " + std::to_string(cellType_53) + "\"" );
    cmd.append(" -n \"cellType_52 " + std::to_string(cellType_52) + "\"" );
    cmd.append(" -n \"cellType_51 " + std::to_string(cellType_51) + "\"" );
    cmd.append(" -n \"cellType_50 " + std::to_string(cellType_50) + "\"" );
    cmd.append(" -n \"cellType_49 " + std::to_string(cellType_49) + "\"" );
    cmd.append(" -n \"cellType_48 " + std::to_string(cellType_48) + "\"" );
    cmd.append(" -n \"cellType_47 " + std::to_string(cellType_47) + "\"" );
    cmd.append(" -n \"cellType_46 " + std::to_string(cellType_46) + "\"" );
    cmd.append(" -n \"cellType_45 " + std::to_string(cellType_45) + "\"" );
    cmd.append(" -n \"cellType_44 " + std::to_string(cellType_44) + "\"" );
    cmd.append(" -n \"cellType_43 " + std::to_string(cellType_43) + "\"" );
    cmd.append(" -n \"cellType_42 " + std::to_string(cellType_42) + "\"" );
    cmd.append(" -n \"cellType_41 " + std::to_string(cellType_41) + "\"" );
    cmd.append(" -n \"cellType_40 " + std::to_string(cellType_40) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );
    
    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_40 % 10) == 0) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 1) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 2) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 3) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 4) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 5) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 6) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 7) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 8) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_40 % 10) == 9) ret += (( nab_40 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_41 % 10) == 0) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 1) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 2) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 3) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 4) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 5) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 6) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 7) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 8) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_41 % 10) == 9) ret += (( nab_41 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_42 % 10) == 0) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 1) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 2) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 3) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 4) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 5) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 6) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 7) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 8) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_42 % 10) == 9) ret += (( nab_42 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_43 % 10) == 0) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 1) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 2) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 3) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 4) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 5) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 6) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 7) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 8) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_43 % 10) == 9) ret += (( nab_43 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_44 % 10) == 0) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 1) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 2) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 3) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 4) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 5) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 6) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 7) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 8) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_44 % 10) == 9) ret += (( nab_44 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_45 % 10) == 0) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 1) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 2) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 3) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 4) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 5) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 6) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 7) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 8) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_45 % 10) == 9) ret += (( nab_45 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_46 % 10) == 0) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 1) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 2) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 3) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 4) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 5) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 6) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 7) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 8) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_46 % 10) == 9) ret += (( nab_46 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_47 % 10) == 0) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 1) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 2) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 3) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 4) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 5) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 6) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 7) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 8) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_47 % 10) == 9) ret += (( nab_47 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_48 % 10) == 0) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 1) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 2) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 3) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 4) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 5) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 6) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 7) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 8) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_48 % 10) == 9) ret += (( nab_48 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_49 % 10) == 0) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 1) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 2) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 3) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 4) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 5) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 6) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 7) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 8) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_49 % 10) == 9) ret += (( nab_49 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_50 % 10) == 0) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 1) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 2) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 3) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 4) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 5) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 6) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 7) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 8) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_50 % 10) == 9) ret += (( nab_50 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_51 % 10) == 0) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 1) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 2) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 3) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 4) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 5) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 6) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 7) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 8) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_51 % 10) == 9) ret += (( nab_51 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_52 % 10) == 0) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 1) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 2) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 3) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 4) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 5) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 6) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 7) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 8) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_52 % 10) == 9) ret += (( nab_52 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_53 % 10) == 0) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 1) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 2) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 3) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 4) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 5) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 6) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 7) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 8) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_53 % 10) == 9) ret += (( nab_53 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_54 % 10) == 0) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 1) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 2) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 3) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 4) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 5) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 6) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 7) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 8) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_54 % 10) == 9) ret += (( nab_54 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

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


    ret /= (19.0 * 16 * FA_TRANSISTOR_COUNT);
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


