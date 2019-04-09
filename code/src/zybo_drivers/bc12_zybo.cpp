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
 * @file   bc12_zybo.cpp
 * @author Andrea Aletto
 * @date   5 apr 2019
 * @brief  Implementation of bc12_zybo class
 ******************************************************************************/

#include "bc12_zybo.h"

static void setComponent(unsigned* bc12_phy_addr, int i, int j, int16_t val);
static int16_t getComponent(unsigned* bc12_phy_addr, int i, int j);

void BC12_zybo::dct(const cv::Mat& input, cv::Mat& output){

    std::cout<<"\nCalling BC12 Zybo Linux driver...";

	int fd;
	unsigned bc12_base_addr = 0x43C00000;
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
	page_addr = (bc12_base_addr & (~(page_size-1)));
	page_offset = bc12_base_addr - page_addr;
	ptr = mmap(NULL, page_size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, page_addr);

	std::cout << "\npage addr: " <<(void*)page_addr;
	std::cout << "\npage offset: " <<(void*)page_offset;
	std::cout << "\nptr: " <<(void*)ptr;
	if (ptr == MAP_FAILED) {
		printf("Mapping indirizzo fisico - indirizzo virtuale FALLITO!\n");
		return;
	}
	
	/* Write values to the device registers */
    for(int i=0; i<input.rows; i++){
        for(int j=0; j<input.cols; j++){
            int16_t val = input.at<int16_t>(i,j);
	        // *((unsigned *)(ptr + (page_offset + 8*i*0x4 + 0x4*j))) = value;		
    		setComponent((unsigned*)(ptr+page_offset),i,j,val);
	        // printf("Going to write onto %08x the value %08x\n", bc12_base_addr, value);
        }
    }

	std::cout<<(void*)(ptr+page_offset + 4*(8*(input.rows-1)+(input.cols-1)));

	// TODO: ci vuole un'attesa?

	/* Read results */
	for(int i=0; i<input.rows; i++){
        for(int j=0; j<input.cols; j++){
            int16_t val = getComponent((unsigned*)(ptr+page_offset),i,j);
			output.at<int16_t>(i,j) = val;
        }
    }

	munmap(ptr, page_size);
}

void setComponent(unsigned* bc12_phy_addr, int i, int j, int16_t val){
	// BC12_AXI_mWriteReg(bc12_phy_addr, 4*(8*i+j), val);
	*((unsigned*)(bc12_phy_addr + 4*(8*i+j))) = (unsigned)val;

}

int16_t getComponent(unsigned* bc12_phy_addr, int i, int j){
	// return (int16_t)(BC12_AXI_mReadReg(bc12_phy_addr, 64*4 + 4*(8*i+j)));
	unsigned* addr = (unsigned*)(bc12_phy_addr + 64*4 + 4*(8*i+j));
	return (int16_t)(*addr);
}