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
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x CB11 -m dssim";
    
    cmd.append(" -n \"nab_77 " + std::to_string(nab_77) + "\"" );
    cmd.append(" -n \"nab_78 " + std::to_string(nab_78) + "\"" );
    cmd.append(" -n \"nab_79 " + std::to_string(nab_79) + "\"" );
    cmd.append(" -n \"nab_80 " + std::to_string(nab_80) + "\"" );
    cmd.append(" -n \"nab_81 " + std::to_string(nab_81) + "\"" );
    cmd.append(" -n \"nab_82 " + std::to_string(nab_82) + "\"" );
    cmd.append(" -n \"nab_83 " + std::to_string(nab_83) + "\"" );
    cmd.append(" -n \"nab_84 " + std::to_string(nab_84) + "\"" );
    cmd.append(" -n \"nab_85 " + std::to_string(nab_85) + "\"" );
    cmd.append(" -n \"nab_86 " + std::to_string(nab_86) + "\"" );
    cmd.append(" -n \"nab_87 " + std::to_string(nab_87) + "\"" );
    cmd.append(" -n \"nab_88 " + std::to_string(nab_88) + "\"" );
    cmd.append(" -n \"nab_89 " + std::to_string(nab_89) + "\"" );
    cmd.append(" -n \"nab_90 " + std::to_string(nab_90) + "\"" );
    cmd.append(" -n \"nab_91 " + std::to_string(nab_91) + "\"" );
    cmd.append(" -n \"nab_92 " + std::to_string(nab_92) + "\"" );
    cmd.append(" -n \"nab_93 " + std::to_string(nab_93) + "\"" );
    cmd.append(" -n \"nab_94 " + std::to_string(nab_94) + "\"" );
    cmd.append(" -n \"nab_95 " + std::to_string(nab_95) + "\"" );
    cmd.append(" -n \"nab_96 " + std::to_string(nab_96) + "\"" );
    cmd.append(" -n \"nab_97 " + std::to_string(nab_97) + "\"" );
    cmd.append(" -n \"nab_98 " + std::to_string(nab_98) + "\"" );
    cmd.append(" -n \"nab_99 " + std::to_string(nab_99) + "\"" );

    cmd.append(" -n \"cellType_77 " + std::to_string(cellType_77) + "\"" );
    cmd.append(" -n \"cellType_78 " + std::to_string(cellType_78) + "\"" );
    cmd.append(" -n \"cellType_79 " + std::to_string(cellType_79) + "\"" );
    cmd.append(" -n \"cellType_80 " + std::to_string(cellType_80) + "\"" );
    cmd.append(" -n \"cellType_81 " + std::to_string(cellType_81) + "\"" );
    cmd.append(" -n \"cellType_82 " + std::to_string(cellType_82) + "\"" );
    cmd.append(" -n \"cellType_83 " + std::to_string(cellType_83) + "\"" );
    cmd.append(" -n \"cellType_84 " + std::to_string(cellType_84) + "\"" );
    cmd.append(" -n \"cellType_85 " + std::to_string(cellType_85) + "\"" );
    cmd.append(" -n \"cellType_86 " + std::to_string(cellType_86) + "\"" );
    cmd.append(" -n \"cellType_87 " + std::to_string(cellType_87) + "\"" );
    cmd.append(" -n \"cellType_88 " + std::to_string(cellType_88) + "\"" );
    cmd.append(" -n \"cellType_89 " + std::to_string(cellType_89) + "\"" );
    cmd.append(" -n \"cellType_90 " + std::to_string(cellType_90) + "\"" );
    cmd.append(" -n \"cellType_91 " + std::to_string(cellType_91) + "\"" );
    cmd.append(" -n \"cellType_92 " + std::to_string(cellType_92) + "\"" );
    cmd.append(" -n \"cellType_93 " + std::to_string(cellType_93) + "\"" );
    cmd.append(" -n \"cellType_94 " + std::to_string(cellType_94) + "\"" );
    cmd.append(" -n \"cellType_95 " + std::to_string(cellType_95) + "\"" );
    cmd.append(" -n \"cellType_96 " + std::to_string(cellType_96) + "\"" );
    cmd.append(" -n \"cellType_97 " + std::to_string(cellType_97) + "\"" );
    cmd.append(" -n \"cellType_98 " + std::to_string(cellType_98) + "\"" );
    cmd.append(" -n \"cellType_99 " + std::to_string(cellType_99) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );    

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_77 % 10) == 0) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 1) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 2) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 3) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 4) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 5) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 6) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 7) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 8) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_77 % 10) == 9) ret += (( nab_77 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_78 % 10) == 0) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 1) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 2) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 3) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 4) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 5) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 6) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 7) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 8) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_78 % 10) == 9) ret += (( nab_78 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_79 % 10) == 0) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 1) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 2) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 3) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 4) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 5) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 6) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 7) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 8) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_79 % 10) == 9) ret += (( nab_79 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_80 % 10) == 0) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 1) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 2) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 3) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 4) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 5) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 6) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 7) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 8) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_80 % 10) == 9) ret += (( nab_80 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_81 % 10) == 0) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 1) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 2) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 3) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 4) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 5) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 6) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 7) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 8) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_81 % 10) == 9) ret += (( nab_81 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_82 % 10) == 0) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 1) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 2) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 3) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 4) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 5) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 6) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 7) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 8) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_82 % 10) == 9) ret += (( nab_82 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_83 % 10) == 0) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 1) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 2) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 3) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 4) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 5) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 6) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 7) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 8) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_83 % 10) == 9) ret += (( nab_83 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_84 % 10) == 0) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 1) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 2) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 3) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 4) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 5) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 6) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 7) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 8) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_84 % 10) == 9) ret += (( nab_84 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_85 % 10) == 0) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 1) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 2) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 3) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 4) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 5) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 6) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 7) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 8) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_85 % 10) == 9) ret += (( nab_85 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_86 % 10) == 0) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 1) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 2) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 3) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 4) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 5) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 6) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 7) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 8) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_86 % 10) == 9) ret += (( nab_86 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_87 % 10) == 0) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 1) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 2) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 3) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 4) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 5) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 6) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 7) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 8) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_87 % 10) == 9) ret += (( nab_87 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_88 % 10) == 0) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 1) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 2) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 3) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 4) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 5) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 6) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 7) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 8) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_88 % 10) == 9) ret += (( nab_88 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_89 % 10) == 0) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 1) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 2) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 3) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 4) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 5) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 6) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 7) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 8) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_89 % 10) == 9) ret += (( nab_89 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_90 % 10) == 0) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 1) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 2) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 3) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 4) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 5) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 6) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 7) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 8) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_90 % 10) == 9) ret += (( nab_90 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_91 % 10) == 0) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 1) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 2) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 3) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 4) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 5) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 6) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 7) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 8) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_91 % 10) == 9) ret += (( nab_91 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_92 % 10) == 0) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 1) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 2) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 3) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 4) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 5) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 6) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 7) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 8) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_92 % 10) == 9) ret += (( nab_92 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_93 % 10) == 0) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 1) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 2) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 3) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 4) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 5) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 6) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 7) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 8) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_93 % 10) == 9) ret += (( nab_93 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_94 % 10) == 0) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 1) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 2) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 3) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 4) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 5) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 6) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 7) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 8) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_94 % 10) == 9) ret += (( nab_94 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_95 % 10) == 0) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 1) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 2) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 3) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 4) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 5) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 6) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 7) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 8) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_95 % 10) == 9) ret += (( nab_95 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_96 % 10) == 0) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 1) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 2) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 3) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 4) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 5) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 6) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 7) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 8) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_96 % 10) == 9) ret += (( nab_96 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_97 % 10) == 0) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 1) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 2) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 3) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 4) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 5) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 6) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 7) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 8) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_97 % 10) == 9) ret += (( nab_97 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    if( (cellType_98 % 10) == 0) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 1) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 2) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 3) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 4) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 5) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 6) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 7) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 8) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_98 % 10) == 9) ret += (( nab_98 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_99 % 10) == 0) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 1) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 2) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 3) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 4) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 5) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 6) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 7) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 8) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_99 % 10) == 9) ret += (( nab_99 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    

    ret /= (23.0 * 16 * FA_TRANSISTOR_COUNT);
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


