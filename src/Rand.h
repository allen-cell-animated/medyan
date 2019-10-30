
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Rand_h
#define MEDYAN_Rand_h

#include <stdio.h>
#include <random>

#include "common.h"

/// A random number generator class.
class Rand {
    
private:
    static uniform_int_distribution<int> _int_distr;
    static uniform_int_distribution<uint64_t> _uint64_t_distr;

public:
    static mt19937 eng;
    #ifdef DEBUGCONSTANTSEED
    static int intcounter;
    static int floatcounter;
    static int chemistrycounter;
	#endif

    ///Get a random floatingpoint between low and high
    static inline floatingpoint randfloatingpoint(floatingpoint low, floatingpoint high) {
        #ifdef DEBUGCONSTANTSEED
        floatcounter++;
		#endif
        return ((float)_int_distr(eng) / numeric_limits<int>::max()) * (high - low) + low;
    }
    ///Get a random integer between low and high
    static inline int randInteger(int low, int high) {
        #ifdef DEBUGCONSTANTSEED
        intcounter++;
        #endif
        int y =_int_distr(eng);
        int x = low + (y % (high - low + 1));
        return x;
    }

    static inline uint64_t randInteger64bit(){
        return _uint64_t_distr(eng);
    }
};

#endif
