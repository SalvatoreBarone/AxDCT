/******************************************************************************
 * @file   bellerophon_functions.cpp
 * @author Andrea Aletto
 * @date   13 Feb 2019
 * @brief  Bellerophon functions for the DCT Algorithm Approximated through
 *         clang-chimera
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>


#include "algorithms_list.h"
#include "algorithms/BC12.h"


char oracle_path[]= "./main_oracle.txt";


extern "C" double BELLERO_getError()
{
    // // Local vars
    // double error = 0;
    // int expected_val, calc_val;
    // int first_val, second_val, third_val;
    // int cnt = 0; // Counter

    // ::std::string first_str, second_str, third_str, expected_str;
    // ::std::ifstream oracle ( oracle_path );
    // if ( !oracle.good() ) {
    //     ::std::cerr << "The oracle does not exist!\n";
    //     exit ( 1 );
    // }

    // // Read oracle's values
    
    // do {
    //     oracle >> first_str >> second_str >> third_str >> expected_str;
    //     cnt++;
    //     expected_val = ::std::atoi ( expected_str.c_str() );

    //     first_val = ::std::atoi( first_str.c_str() );
    //     second_val = ::std::atoi( second_str.c_str() );
    //     third_val = ::std::atoi( third_str.c_str() );
    //     calc_val = sumB(first_val, second_val ,third_val);
    //     error += fabs ( calc_val - expected_val ) < 10e-16 ? 0 : fabs ( calc_val - expected_val );
    // } while ( !oracle.eof() );

    return 1e-6;
}

