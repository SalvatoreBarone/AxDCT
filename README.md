# AxDCT
[![Build Status](https://travis-ci.com/andreaaletto/AxDCT.svg?token=vSvxrpZbWB2t5qeWZ4CJ&branch=master)](https://travis-ci.com/andreaaletto/AxDCT) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

AxDCT is a collection of inexact DCT algorithms useful to evaluate the effect of _Approximate Computing_ techniques on _Image Computing_.
This software suite is meant be used together with [IIDEAA](http://wpage.unina.it/mario.barbareschi/old/iideaa/handson/) tools (_Clang-Chimera_ and _Bellerophon_), but, however, it is possible to use it standalone.

## Build and install

In order to build **AxDCT** source code you will need an installation of a Linux distro (in the following Ubuntu 18.04 will be assumed) with the following dependencies:
- An up-to-date development software suite (_build-essential, git, cmake, make, etc._),
- An installation of [OpenCV](https://opencv.org/) (minimum supported version 4.0.1),
- An installation of [CNL](https://github.com/johnmcfarlane/cnl) Library for fixed point math.

In order to build and run code mutation stuff, you will need the following dependencies, too:
- An installation of LLVM/Clang 3.9.1 as explaind in [clang-chimera](https://github.com/andreaaletto/clang-chimera) repository.
- An installation of clang-chimera
- An installation of [Bellerophon](https://github.com/andreaaletto/Bellerophon)
- An installation of [AxC-Adders](https://github.com/andreaaletto/AxC-Adders) static library

### Dependencies installation
Firstly, you need to update your development machine. To do so, run the following commands:
```
# apt-get update
# apt-get upgrade -y
# apt-get install git make ninja-build build-essential wget unzip zlib1g-dev libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev python-dev python-numpy libtbb2 libtbb-dev libjpeg-dev libpng-dev libtiff-dev jasper libdc1394-22-dev libffi-dev libedit-dev libncurses5-dev libboost-dev
```

In order to configure AxDCT project it is required the latest version of cmake (3.14), so you need to remove it (if already installed) and to install the newer version:
```
# apt remove --purge --auto-remove cmake
# apt purge --auto-remove cmake
$ wget https://github.com/Kitware/CMake/releases/download/v3.14.0/cmake-3.14.0.tar.gz
$ tar xvfz cmake-3.14.0.tar.gz
$ cd cmake-3.14.0
$ ./bootstrap
$ make -j4 && sudo make install
$ cmake --version
```

Now download, build and install OpenCV library:
```
$ wget https://github.com/opencv/opencv/archive/4.0.1.zip
$ unzip 4.0.1.zip && mv opencv-4.0.1 opencv
$ cd opencv && mkdir build && cd build
$ cmake -D CMAKE_BUILD_TYPE=RELEASE -D BUILD_EXAMPLES=OFF -D BUILD_SHARED_LIBS=OFF -D BUILD_opencv_apps=OFF -D BUILD_DOCS=OFF -D BUILD_PERF_TESTS=OFF -D BUILD_TESTS=OFF -D CMAKE_INSTALL_PREFIX=/usr/local ..
$ make -j4 
# make install
# sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/opencv.conf'
# ldconfig
```

Now download and install CNL library:
```
$ git clone https://github.com/johnmcfarlane/cnl.git && cd cnl
$ mkdir build && cd build
$ cmake .. -G Ninja
# ninja install
```

### Additional dependencies installation
If you want to evaluate AxC techniques' impact onto AxDCT algorithms, you need to install IIDEAA tools. Please refer to [clang-chimera](https://github.com/andreaaletto/clang-chimera) readme to discover how to correctly install LLVM/Clang 3.9.1 and clang-chimera tool, then install [Bellerophon](https://github.com/andreaaletto/Bellerophon) following the readme. At last, install [AxC-Adders](https://github.com/andreaaletto/AxC-Adders) library.
To test if all the additional dependencies are been correctly installed, do the following commands: 
```
$ llvm-config --version
$ clang --version
$ clang-chimera -version
$ bellerophon -version
$ ls /usr/lib/libInexactAdders.a
```

### Compile AxDCT
--------

With all dependencies installed, now you can download, build and run AxDCT:
```
$ git clone https://github.com/andreaaletto/AxDCT.git && cd AxDCT/code
$ mkdir build && cd build
$ cmake .. -G Ninja
$ ninja
```

At the end of the process, you will find the following executables in ```bin``` folder:
* ```axdct``` - To obtain and visualize an image direct-transformed with an approximate DCT and inverse-transformed with an exact IDCT.
* ```image_eval``` - To calculate a particular metric upon an image direct-transformed with an approximate DCT and inverse-transformed with an exact IDCT.

Running ```axdct -h``` or ```image_eval -h```, it is possible to discover all the option flag to pass to the program.

#### A simple example
Just for example purposes, after compilation, we run the -l option order to discover the supported algorithms:
```
$ bin/axdct -l

Output:
Supported algorithms:

- BC12		
- CB11		
- BAS08	
- BAS09	
- BAS11 (a=0.0)
- BAS11 (a=0.5)
- BAS11 (a=1.0)
- BAS11 (a=2.0)
- PEA12	
- PEA14	

Supported metrics:

- PSNR
- MSE
- MD
- AD
- MSSIM
```

Now we decide to see the effect of the algorithm BC12, so we run:
```
$ bin/axdct -i /path/to/image.bmp -x bc12 
```

Using OpenCV library the following output is shown:
TODO: screenshot

If now we want to evaluate the PSNR of that process, we run:
```
$ bin/image_eval -i /path/to/image.bmp -x bc12 -m PSNR

Output:
************** PSNR **************

    BC12         	26.015

**********************************
```

### Prepare code-mutants with clang-chimera
If you want to evaluate code mutation on AxDCT, you will need to generate mutants. So, assuming you've already gone through _Additional dependencies installation_ section, run the following commands:
```
$ cd $HOME/AxDCT/chimera
$ chmod +x launch_chimera.sh
$ ./launch_chimera.sh
```

This will generate in the folder ```output``` a new .cpp file for each algorithm implementation, substituting every sum operation in the AxDCT algorithm, with a function call to the model of inexact-addition of ```AxC-Adders``` library.

### Compile AxDCT with mutants
In order to include mutants in the build process, you need to re-run cmake as follows:
```
$ cd $HOME/AxDCT/code/build
$ cmake .. -G Ninja -DAXDCT_BUILD_MUTANTS=ON
$ ninja
```

This compile process will generate a new executable, too:
* ```mutants_eval``` - To evaluate metrics onto the AxDCT mutation.

In particular, this executable introduces the flag -n, used to assign a static value to a global variable. The list of the available global variables can be found inside the csv produced by clang-chimera.

### Compile ZyBo Drivers
It is possible to build drivers for Linux as MMAP on /dev/mem, meant to be used together with the hardware design of [AxC-Adders_vhdl](https://github.com/andreaaletto/AxC-Adders_vhdl), implemented as AXI peripheral on armhf architecture.

In order to include zybo drivers in the build process, you need to re-run cmake as follows:
```
$ cd $HOME/AxDCT/code/build
$ cmake .. -G Ninja -DAXDCT_BUILD_ZYBO=ON
$ ninja
```

TODO: documentation

#### LICENSE
--------

* [GPLV3.0](https://www.gnu.org/licenses/licenses.html)

#### Contributing
----------

Github is for social coding: if you want to write code, I encourage contributions through pull requests from forks of this repository. 
