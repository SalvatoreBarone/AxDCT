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
 * @file   pea14_zybo.cpp
 * @author Andrea Aletto
 * @date   5 apr 2019
 * @brief  Implementation of pea14_zybo class
 ******************************************************************************/

#include "pea14_zybo.h"

//#define __DEBUG__

void PEA14_zybo::dct(const cv::Mat& input, cv::Mat& output){
#ifdef __DEBUG__
    std::cout<<"\nCalling PEA14 Zybo Linux driver...";
#endif

	int fd;
	unsigned pea14_base_addr = 0x43C00000; //BE SURE THAT THIS IS THE BASE ADDRESS OF THE AXI PERIPHERAL IN THE HW DESIGN
	unsigned page_addr, page_offset;
	void *ptr;
	unsigned page_size=sysconf(_SC_PAGESIZE);

    /* Open /dev/mem file */
	fd = open ("/dev/mem", O_RDWR);
	if (fd < 1) {
		std::cerr << "\nERROR: Cannot open /dev/mem";
		return;
	}

    /* mmap the device into memory */
	page_addr = (pea14_base_addr & (~(page_size-1)));
	page_offset = pea14_base_addr - page_addr;
	ptr = mmap(NULL, page_size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, page_addr);

#ifdef __DEBUG__
	std::cout << "\npage addr: " <<(void*)page_addr;
	std::cout << "\npage offset: " <<(void*)page_offset;
	std::cout << "\npage size: " <<(void*)page_size;
	std::cout << "\nptr: " <<(void*)ptr;
	std::cout << "\nslv_reg63 addr: " << (void*)(ptr+page_offset + 4*(8*(input.rows-1)+(input.cols-1)));
	std::cout << "\nslv_reg127 addr: " << (void*)(ptr+page_offset + 64*4 + 4*(8*(input.rows-1)+(input.cols-1)));
#endif

	if (ptr == MAP_FAILED) {
		printf("Mapping indirizzo fisico - indirizzo virtuale FALLITO!\n");
		return;
	}
	
	/* Write values to the device registers */
    for(int i=0; i<input.rows; i++){
        for(int j=0; j<input.cols; j++){
            int16_t val = input.at<int16_t>(i,j);
			*((uint32_t*)(ptr+page_offset + 4*(8*i+j))) = (uint32_t)val;
        }
    }

	/* Read results */
	for(int i=0; i<input.rows; i++){
        for(int j=0; j<input.cols; j++){
			uint32_t* addr = (uint32_t*)(ptr+page_offset + 64*4 + 4*(8*i+j));
			output.at<int16_t>(i,j) = (int16_t)(*addr);
        }
    }

	munmap(ptr, page_size);
	close(fd);
}
