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
    std::string cmd = "~/AxDCT/code/build/bin/mutant_eval -s -f ~/image_dataset -x PEA12 -m dssim";
    
    cmd.append( " -n \"nab_100 " + std::to_string(nab_100) + "\"" );
    cmd.append( " -n \"nab_101 " + std::to_string(nab_101) + "\"" );
    cmd.append( " -n \"nab_102 " + std::to_string(nab_102) + "\"" );
    cmd.append( " -n \"nab_103 " + std::to_string(nab_103) + "\"" );
    cmd.append( " -n \"nab_104 " + std::to_string(nab_104) + "\"" );
    cmd.append( " -n \"nab_105 " + std::to_string(nab_105) + "\"" );
    cmd.append( " -n \"nab_106 " + std::to_string(nab_106) + "\"" );
    cmd.append( " -n \"nab_107 " + std::to_string(nab_107) + "\"" );
    cmd.append( " -n \"nab_108 " + std::to_string(nab_108) + "\"" );
    cmd.append( " -n \"nab_109 " + std::to_string(nab_109) + "\"" );
    cmd.append( " -n \"nab_110 " + std::to_string(nab_110) + "\"" );
    cmd.append( " -n \"nab_111 " + std::to_string(nab_111) + "\"" );
    cmd.append( " -n \"nab_112 " + std::to_string(nab_112) + "\"" );
    cmd.append( " -n \"nab_113 " + std::to_string(nab_113) + "\"" );
    cmd.append( " -n \"nab_114 " + std::to_string(nab_114) + "\"" );
    cmd.append( " -n \"nab_115 " + std::to_string(nab_115) + "\"" );
    cmd.append( " -n \"nab_116 " + std::to_string(nab_116) + "\"" );
    cmd.append( " -n \"nab_117 " + std::to_string(nab_117) + "\"" );
    cmd.append( " -n \"nab_118 " + std::to_string(nab_118) + "\"" );
    cmd.append( " -n \"nab_119 " + std::to_string(nab_119) + "\"" );
    cmd.append( " -n \"nab_120 " + std::to_string(nab_120) + "\"" );
    cmd.append( " -n \"nab_121 " + std::to_string(nab_121) + "\"" );
    cmd.append( " -n \"nab_122 " + std::to_string(nab_122) + "\"" );
    cmd.append( " -n \"nab_123 " + std::to_string(nab_123) + "\"" );

    cmd.append(" -n \"cellType_100 " + std::to_string(cellType_100) + "\"" );
    cmd.append(" -n \"cellType_101 " + std::to_string(cellType_101) + "\"" );
    cmd.append(" -n \"cellType_102 " + std::to_string(cellType_102) + "\"" );
    cmd.append(" -n \"cellType_103 " + std::to_string(cellType_103) + "\"" );
    cmd.append(" -n \"cellType_104 " + std::to_string(cellType_104) + "\"" );
    cmd.append(" -n \"cellType_105 " + std::to_string(cellType_105) + "\"" );
    cmd.append(" -n \"cellType_106 " + std::to_string(cellType_106) + "\"" );
    cmd.append(" -n \"cellType_107 " + std::to_string(cellType_107) + "\"" );
    cmd.append(" -n \"cellType_108 " + std::to_string(cellType_108) + "\"" );
    cmd.append(" -n \"cellType_109 " + std::to_string(cellType_109) + "\"" );
    cmd.append(" -n \"cellType_110 " + std::to_string(cellType_110) + "\"" );
    cmd.append(" -n \"cellType_111 " + std::to_string(cellType_111) + "\"" );
    cmd.append(" -n \"cellType_112 " + std::to_string(cellType_112) + "\"" );
    cmd.append(" -n \"cellType_113 " + std::to_string(cellType_113) + "\"" );
    cmd.append(" -n \"cellType_114 " + std::to_string(cellType_114) + "\"" );
    cmd.append(" -n \"cellType_115 " + std::to_string(cellType_115) + "\"" );
    cmd.append(" -n \"cellType_116 " + std::to_string(cellType_116) + "\"" );
    cmd.append(" -n \"cellType_117 " + std::to_string(cellType_117) + "\"" );
    cmd.append(" -n \"cellType_118 " + std::to_string(cellType_118) + "\"" );
    cmd.append(" -n \"cellType_119 " + std::to_string(cellType_119) + "\"" );
    cmd.append(" -n \"cellType_120 " + std::to_string(cellType_120) + "\"" );
    cmd.append(" -n \"cellType_121 " + std::to_string(cellType_121) + "\"" );
    cmd.append(" -n \"cellType_122 " + std::to_string(cellType_122) + "\"" );
    cmd.append(" -n \"cellType_123 " + std::to_string(cellType_123) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );
    

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_100 % 10) == 0) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 1) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 2) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 3) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 4) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 5) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 6) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 7) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 8) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_100 % 10) == 9) ret += (( nab_100 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_101 % 10) == 0) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 1) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 2) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 3) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 4) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 5) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 6) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 7) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 8) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_101 % 10) == 9) ret += (( nab_101 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_102 % 10) == 0) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 1) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 2) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 3) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 4) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 5) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 6) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 7) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 8) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_102 % 10) == 9) ret += (( nab_102 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_103 % 10) == 0) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 1) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 2) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 3) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 4) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 5) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 6) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 7) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 8) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_103 % 10) == 9) ret += (( nab_103 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_104 % 10) == 0) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 1) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 2) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 3) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 4) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 5) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 6) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 7) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 8) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_104 % 10) == 9) ret += (( nab_104 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_105 % 10) == 0) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 1) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 2) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 3) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 4) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 5) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 6) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 7) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 8) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_105 % 10) == 9) ret += (( nab_105 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_106 % 10) == 0) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 1) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 2) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 3) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 4) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 5) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 6) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 7) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 8) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_106 % 10) == 9) ret += (( nab_106 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_107 % 10) == 0) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 1) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 2) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 3) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 4) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 5) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 6) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 7) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 8) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_107 % 10) == 9) ret += (( nab_107 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_108 % 10) == 0) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 1) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 2) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 3) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 4) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 5) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 6) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 7) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 8) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_108 % 10) == 9) ret += (( nab_108 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_109 % 10) == 0) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 1) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 2) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 3) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 4) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 5) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 6) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 7) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 8) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_109 % 10) == 9) ret += (( nab_109 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_110 % 10) == 0) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 1) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 2) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 3) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 4) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 5) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 6) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 7) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 8) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_110 % 10) == 9) ret += (( nab_110 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_111 % 10) == 0) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 1) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 2) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 3) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 4) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 5) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 6) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 7) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 8) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_111 % 10) == 9) ret += (( nab_111 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_112 % 10) == 0) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 1) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 2) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 3) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 4) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 5) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 6) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 7) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 8) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_112 % 10) == 9) ret += (( nab_112 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_113 % 10) == 0) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 1) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 2) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 3) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 4) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 5) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 6) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 7) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 8) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_113 % 10) == 9) ret += (( nab_113 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_114 % 10) == 0) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 1) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 2) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 3) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 4) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 5) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 6) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 7) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 8) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_114 % 10) == 9) ret += (( nab_114 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_115 % 10) == 0) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 1) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 2) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 3) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 4) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 5) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 6) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 7) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 8) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_115 % 10) == 9) ret += (( nab_115 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_116 % 10) == 0) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 1) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 2) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 3) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 4) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 5) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 6) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 7) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 8) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_116 % 10) == 9) ret += (( nab_116 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_117 % 10) == 0) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 1) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 2) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 3) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 4) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 5) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 6) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 7) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 8) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_117 % 10) == 9) ret += (( nab_117 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_118 % 10) == 0) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 1) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 2) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 3) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 4) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 5) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 6) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 7) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 8) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_118 % 10) == 9) ret += (( nab_118 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_119 % 10) == 0) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 1) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 2) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 3) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 4) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 5) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 6) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 7) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 8) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_119 % 10) == 9) ret += (( nab_119 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_120 % 10) == 0) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 1) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 2) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 3) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 4) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 5) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 6) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 7) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 8) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_120 % 10) == 9) ret += (( nab_120 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_121 % 10) == 0) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 1) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 2) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 3) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 4) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 5) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 6) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 7) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 8) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_121 % 10) == 9) ret += (( nab_121 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_122 % 10) == 0) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 1) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 2) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 3) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 4) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 5) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 6) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 7) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 8) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_122 % 10) == 9) ret += (( nab_122 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_123 % 10) == 0) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 1) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 2) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 3) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 4) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 5) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 6) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 7) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 8) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_123 % 10) == 9) ret += (( nab_123 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );
    
    ret /= (24.0 * 16 * FA_TRANSISTOR_COUNT);
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


