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
    
    cmd.append(" -n \"nab_124 " + std::to_string(nab_124) + "\"" );
    cmd.append(" -n \"nab_125 " + std::to_string(nab_125) + "\"" );
    cmd.append(" -n \"nab_126 " + std::to_string(nab_126) + "\"" );
    cmd.append(" -n \"nab_127 " + std::to_string(nab_127) + "\"" );
    cmd.append(" -n \"nab_128 " + std::to_string(nab_128) + "\"" );
    cmd.append(" -n \"nab_129 " + std::to_string(nab_129) + "\"" );
    cmd.append(" -n \"nab_130 " + std::to_string(nab_130) + "\"" );
    cmd.append(" -n \"nab_131 " + std::to_string(nab_131) + "\"" );
    cmd.append(" -n \"nab_132 " + std::to_string(nab_132) + "\"" );
    cmd.append(" -n \"nab_133 " + std::to_string(nab_133) + "\"" );
    cmd.append(" -n \"nab_134 " + std::to_string(nab_134) + "\"" );
    cmd.append(" -n \"nab_135 " + std::to_string(nab_135) + "\"" );
    cmd.append(" -n \"nab_136 " + std::to_string(nab_136) + "\"" );
    cmd.append(" -n \"nab_137 " + std::to_string(nab_137) + "\"" );
    cmd.append(" -n \"nab_138 " + std::to_string(nab_138) + "\"" );


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
    cmd.append(" -n \"cellType_69 " + std::to_string(cellType_69) + "\"" );
    cmd.append(" -n \"cellType_70 " + std::to_string(cellType_70) + "\"" );
    cmd.append(" -n \"cellType_71 " + std::to_string(cellType_71) + "\"" );
    cmd.append(" -n \"cellType_72 " + std::to_string(cellType_72) + "\"" );
    cmd.append(" -n \"cellType_73 " + std::to_string(cellType_73) + "\"" );
    cmd.append(" -n \"cellType_74 " + std::to_string(cellType_74) + "\"" );
    cmd.append(" -n \"cellType_75 " + std::to_string(cellType_75) + "\"" );
    cmd.append(" -n \"cellType_76 " + std::to_string(cellType_76) + "\"" );

    cmd.append(" -n \"base_0 " + std::to_string(base_0) + "\"" );

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
}

extern "C" double BELLERO_Reward()
{
    double ret = 0.0;

    if( (cellType_124 % 10) == 0) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 1) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 2) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 3) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 4) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 5) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 6) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 7) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 8) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_124 % 10) == 9) ret += (( nab_124 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_125 % 10) == 0) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 1) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 2) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 3) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 4) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 5) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 6) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 7) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 8) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_125 % 10) == 9) ret += (( nab_125 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_126 % 10) == 0) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 1) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 2) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 3) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 4) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 5) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 6) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 7) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 8) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_126 % 10) == 9) ret += (( nab_126 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_127 % 10) == 0) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 1) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 2) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 3) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 4) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 5) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 6) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 7) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 8) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_127 % 10) == 9) ret += (( nab_127 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_128 % 10) == 0) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 1) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 2) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 3) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 4) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 5) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 6) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 7) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 8) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_128 % 10) == 9) ret += (( nab_128 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_129 % 10) == 0) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 1) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 2) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 3) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 4) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 5) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 6) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 7) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 8) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_129 % 10) == 9) ret += (( nab_129 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_130 % 10) == 0) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 1) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 2) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 3) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 4) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 5) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 6) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 7) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 8) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_130 % 10) == 9) ret += (( nab_130 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_131 % 10) == 0) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 1) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 2) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 3) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 4) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 5) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 6) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 7) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 8) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_131 % 10) == 9) ret += (( nab_131 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_132 % 10) == 0) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 1) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 2) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 3) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 4) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 5) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 6) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 7) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 8) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_132 % 10) == 9) ret += (( nab_132 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_133 % 10) == 0) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 1) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 2) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 3) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 4) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 5) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 6) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 7) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 8) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_133 % 10) == 9) ret += (( nab_133 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_134 % 10) == 0) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 1) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 2) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 3) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 4) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 5) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 6) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 7) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 8) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_134 % 10) == 9) ret += (( nab_134 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_135 % 10) == 0) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 1) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 2) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 3) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 4) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 5) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 6) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 7) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 8) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_135 % 10) == 9) ret += (( nab_135 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_136 % 10) == 0) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 1) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 2) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 3) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 4) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 5) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 6) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 7) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 8) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_136 % 10) == 9) ret += (( nab_136 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_137 % 10) == 0) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 1) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 2) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 3) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 4) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 5) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 6) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 7) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 8) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_137 % 10) == 9) ret += (( nab_137 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );

    if( (cellType_138 % 10) == 0) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - InAx1_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 1) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - InAx2_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 2) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - InAx3_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 3) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AMA1_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 4) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AMA2_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 5) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AMA3_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 6) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AMA4_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 7) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AXA1_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 8) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AXA2_TRANSISTOR_COUNT) ) );
    if( (cellType_138 % 10) == 9) ret += (( nab_138 * (FA_TRANSISTOR_COUNT - AXA3_TRANSISTOR_COUNT) ) );


    ret /= (15.0 * 16 * FA_TRANSISTOR_COUNT);
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


