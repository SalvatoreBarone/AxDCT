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

#include "bellero_nablist.h"

std::string exec(const char* cmd);

extern "C" double BELLERO_getError()
{
    std::string cmd = "/home/andrea/AxDCT/code/build/bin/psnr_eval -s -i /home/andrea/lena.bmp -x BC12";
    
    // cmd.append(" -n \"nab_0 " + std::to_string(nab_0) + "\"" );
    // cmd.append(" -n \"nab_1 " + std::to_string(nab_1) + "\"" );
    // cmd.append(" -n \"nab_2 " + std::to_string(nab_2) + "\"" );
    // cmd.append(" -n \"nab_3 " + std::to_string(nab_3) + "\"" );
    // cmd.append(" -n \"nab_4 " + std::to_string(nab_4) + "\"" );
    // cmd.append(" -n \"nab_5 " + std::to_string(nab_5) + "\"" );
    // cmd.append(" -n \"nab_6 " + std::to_string(nab_6) + "\"" );
    // cmd.append(" -n \"nab_7 " + std::to_string(nab_7) + "\"" );
    // cmd.append(" -n \"nab_8 " + std::to_string(nab_8) + "\"" );
    // cmd.append(" -n \"nab_9 " + std::to_string(nab_9) + "\"" );
    // cmd.append(" -n \"nab_10 " + std::to_string(nab_10) + "\"" );
    // cmd.append(" -n \"nab_11 " + std::to_string(nab_11) + "\"" );
    // cmd.append(" -n \"nab_12 " + std::to_string(nab_12) + "\"" );
    // cmd.append(" -n \"nab_13 " + std::to_string(nab_13) + "\"" );
    // cmd.append(" -n \"nab_14 " + std::to_string(nab_14) + "\"" );
    // cmd.append(" -n \"nab_15 " + std::to_string(nab_15) + "\"" );
    // cmd.append(" -n \"nab_16 " + std::to_string(nab_16) + "\"" );
    // cmd.append(" -n \"nab_17 " + std::to_string(nab_17) + "\"" );
    // cmd.append(" -n \"nab_18 " + std::to_string(nab_18) + "\"" );
    // cmd.append(" -n \"nab_19 " + std::to_string(nab_19) + "\"" );
    // cmd.append(" -n \"nab_20 " + std::to_string(nab_20) + "\"" );
    // cmd.append(" -n \"nab_21 " + std::to_string(nab_21) + "\"" );
    // cmd.append(" -n \"nab_22 " + std::to_string(nab_22) + "\"" );
    // cmd.append(" -n \"nab_23 " + std::to_string(nab_23) + "\"" );
    // cmd.append(" -n \"nab_24 " + std::to_string(nab_24) + "\"" );
    // cmd.append(" -n \"nab_25 " + std::to_string(nab_25) + "\"" );
    // cmd.append(" -n \"nab_26 " + std::to_string(nab_26) + "\"" );
    // cmd.append(" -n \"nab_27 " + std::to_string(nab_27) + "\"" );
    // cmd.append(" -n \"nab_28 " + std::to_string(nab_28) + "\"" );
    // cmd.append(" -n \"nab_29 " + std::to_string(nab_29) + "\"" );
    // cmd.append(" -n \"nab_30 " + std::to_string(nab_30) + "\"" );
    // cmd.append(" -n \"nab_31 " + std::to_string(nab_31) + "\"" );
    // cmd.append(" -n \"nab_32 " + std::to_string(nab_32) + "\"" );
    // cmd.append(" -n \"nab_33 " + std::to_string(nab_33) + "\"" );
    // cmd.append(" -n \"nab_34 " + std::to_string(nab_34) + "\"" );
    // cmd.append(" -n \"nab_35 " + std::to_string(nab_35) + "\"" );
    // cmd.append(" -n \"nab_36 " + std::to_string(nab_36) + "\"" );
    // cmd.append(" -n \"nab_37 " + std::to_string(nab_37) + "\"" );
    // cmd.append(" -n \"nab_38 " + std::to_string(nab_38) + "\"" );
    // cmd.append(" -n \"nab_39 " + std::to_string(nab_39) + "\"" );
    // cmd.append(" -n \"nab_40 " + std::to_string(nab_40) + "\"" );
    // cmd.append(" -n \"nab_41 " + std::to_string(nab_41) + "\"" );
    // cmd.append(" -n \"nab_42 " + std::to_string(nab_42) + "\"" );
    // cmd.append(" -n \"nab_43 " + std::to_string(nab_43) + "\"" );
    // cmd.append(" -n \"nab_44 " + std::to_string(nab_44) + "\"" );
    // cmd.append(" -n \"nab_45 " + std::to_string(nab_45) + "\"" );
    // cmd.append(" -n \"nab_46 " + std::to_string(nab_46) + "\"" );
    // cmd.append(" -n \"nab_47 " + std::to_string(nab_47) + "\"" );
    // cmd.append(" -n \"nab_48 " + std::to_string(nab_48) + "\"" );
    // cmd.append(" -n \"nab_49 " + std::to_string(nab_49) + "\"" );
    // cmd.append(" -n \"nab_50 " + std::to_string(nab_50) + "\"" );
    // cmd.append(" -n \"nab_51 " + std::to_string(nab_51) + "\"" );
    // cmd.append(" -n \"nab_52 " + std::to_string(nab_52) + "\"" );
    // cmd.append(" -n \"nab_53 " + std::to_string(nab_53) + "\"" );
    // cmd.append(" -n \"nab_54 " + std::to_string(nab_54) + "\"" );
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
    // cmd.append(" -n \"nab_69 " + std::to_string(nab_69) + "\"" );
    // cmd.append(" -n \"nab_70 " + std::to_string(nab_70) + "\"" );
    // cmd.append(" -n \"nab_71 " + std::to_string(nab_71) + "\"" );
    // cmd.append(" -n \"nab_72 " + std::to_string(nab_72) + "\"" );
    // cmd.append(" -n \"nab_73 " + std::to_string(nab_73) + "\"" );
    // cmd.append(" -n \"nab_74 " + std::to_string(nab_74) + "\"" );
    // cmd.append(" -n \"nab_75 " + std::to_string(nab_75) + "\"" );
    // cmd.append(" -n \"nab_76 " + std::to_string(nab_76) + "\"" );
    // cmd.append(" -n \"nab_77 " + std::to_string(nab_77) + "\"" );
    // cmd.append(" -n \"nab_78 " + std::to_string(nab_78) + "\"" );
    // cmd.append(" -n \"nab_79 " + std::to_string(nab_79) + "\"" );
    // cmd.append(" -n \"nab_80 " + std::to_string(nab_80) + "\"" );
    // cmd.append(" -n \"nab_81 " + std::to_string(nab_81) + "\"" );
    // cmd.append(" -n \"nab_82 " + std::to_string(nab_82) + "\"" );
    // cmd.append(" -n \"nab_83 " + std::to_string(nab_83) + "\"" );
    // cmd.append(" -n \"nab_84 " + std::to_string(nab_84) + "\"" );
    // cmd.append(" -n \"nab_85 " + std::to_string(nab_85) + "\"" );
    // cmd.append(" -n \"nab_86 " + std::to_string(nab_86) + "\"" );
    // cmd.append(" -n \"nab_87 " + std::to_string(nab_87) + "\"" );
    // cmd.append(" -n \"nab_88 " + std::to_string(nab_88) + "\"" );
    // cmd.append(" -n \"nab_89 " + std::to_string(nab_89) + "\"" );
    // cmd.append(" -n \"nab_90 " + std::to_string(nab_90) + "\"" );
    // cmd.append(" -n \"nab_91 " + std::to_string(nab_91) + "\"" );
    // cmd.append(" -n \"nab_92 " + std::to_string(nab_92) + "\"" );
    // cmd.append(" -n \"nab_93 " + std::to_string(nab_93) + "\"" );
    // cmd.append(" -n \"nab_94 " + std::to_string(nab_94) + "\"" );
    // cmd.append(" -n \"nab_95 " + std::to_string(nab_95) + "\"" );
    // cmd.append(" -n \"nab_96 " + std::to_string(nab_96) + "\"" );
    // cmd.append(" -n \"nab_97 " + std::to_string(nab_97) + "\"" );
    // cmd.append(" -n \"nab_98 " + std::to_string(nab_98) + "\"" );
    // cmd.append(" -n \"nab_99 " + std::to_string(nab_99) + "\"" );
    // cmd.append(" -n \"nab_100 " + std::to_string(nab_100) + "\"" );
    // cmd.append(" -n \"nab_101 " + std::to_string(nab_101) + "\"" );
    // cmd.append(" -n \"nab_102 " + std::to_string(nab_102) + "\"" );
    // cmd.append(" -n \"nab_103 " + std::to_string(nab_103) + "\"" );
    // cmd.append(" -n \"nab_104 " + std::to_string(nab_104) + "\"" );
    // cmd.append(" -n \"nab_105 " + std::to_string(nab_105) + "\"" );
    // cmd.append(" -n \"nab_106 " + std::to_string(nab_106) + "\"" );
    // cmd.append(" -n \"nab_107 " + std::to_string(nab_107) + "\"" );
    // cmd.append(" -n \"nab_108 " + std::to_string(nab_108) + "\"" );
    // cmd.append(" -n \"nab_109 " + std::to_string(nab_109) + "\"" );
    // cmd.append(" -n \"nab_110 " + std::to_string(nab_110) + "\"" );
    // cmd.append(" -n \"nab_111 " + std::to_string(nab_111) + "\"" );
    // cmd.append(" -n \"nab_112 " + std::to_string(nab_112) + "\"" );
    // cmd.append(" -n \"nab_113 " + std::to_string(nab_113) + "\"" );
    // cmd.append(" -n \"nab_114 " + std::to_string(nab_114) + "\"" );
    // cmd.append(" -n \"nab_115 " + std::to_string(nab_115) + "\"" );
    // cmd.append(" -n \"nab_116 " + std::to_string(nab_116) + "\"" );
    // cmd.append(" -n \"nab_117 " + std::to_string(nab_117) + "\"" );
    // cmd.append(" -n \"nab_118 " + std::to_string(nab_118) + "\"" );
    // cmd.append(" -n \"nab_119 " + std::to_string(nab_119) + "\"" );
    // cmd.append(" -n \"nab_120 " + std::to_string(nab_120) + "\"" );
    // cmd.append(" -n \"nab_121 " + std::to_string(nab_121) + "\"" );
    // cmd.append(" -n \"nab_122 " + std::to_string(nab_122) + "\"" );
    // cmd.append(" -n \"nab_123 " + std::to_string(nab_123) + "\"" );
    // cmd.append(" -n \"nab_124 " + std::to_string(nab_124) + "\"" );
    // cmd.append(" -n \"nab_125 " + std::to_string(nab_125) + "\"" );
    // cmd.append(" -n \"nab_126 " + std::to_string(nab_126) + "\"" );
    // cmd.append(" -n \"nab_127 " + std::to_string(nab_127) + "\"" );
    // cmd.append(" -n \"nab_128 " + std::to_string(nab_128) + "\"" );

    std::string retstring = exec(cmd.c_str());
    return stod(retstring);
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


