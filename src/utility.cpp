
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "utility.h"

#if defined(_MSC_VER)
// preprocessors to use __rdtsc() for MSVC only
#  include <intrin.h>
#  pragma intrinsic(__rdtsc)
#endif

unsigned long long rdtsc(){

#if defined(_MSC_VER) // MSVC only
    return __rdtsc();

#else
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;

#endif
}
