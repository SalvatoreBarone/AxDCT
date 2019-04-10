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
 * @file   bas09_zybo.h
 * @author Andrea Aletto
 * @date   10 apr 2019
 * @brief  Declaration of bas09_zybo class
 ******************************************************************************/

#ifndef _BAS09_ZYBO_DRIVER_H
#define _BAS09_ZYBO_DRIVER_H

#include <opencv2/opencv.hpp>
#include "AxDCT_algorithm.h"
#include "BAS09.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>

class BAS09_zybo : public BAS09
{

    private:
        void dct(const cv::Mat& input, cv::Mat& output);

    public:
        BAS09_zybo() : BAS09() {};
};

#endif /* _BAS09_ZYBO_DRIVER_H */