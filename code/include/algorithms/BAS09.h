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
 * @file   BAS09.h
 * @author Andrea Aletto
 * @date   11 feb 2019
 * @brief  Declaration of BAS09 algorithm class
 ******************************************************************************/

#ifndef _BAS09_H
#define _BAS09_H

#include <opencv2/opencv.hpp>

    class BAS09
    {
        private:
            BAS09(){};

        public:
            static void retrieveParameters(cv::Mat& T, cv::Mat& D, cv::Mat& Q, cv::Mat& CQ);
            static cv::Mat getT();
            static cv::Mat getD();
            static cv::Mat getQ();
            static cv::Mat getCQ();
    };

#endif /* _BAS09_H */