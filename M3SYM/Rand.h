
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Rand_h
#define M3SYM_Rand_h

#include <stdio.h>
#include <random>

#include "common.h"

/// A random number generator class.
class Rand {
    
private:
    static mt19937 _eng;
    static uniform_int_distribution<int> _int_distr;
    
public:
    ///Get a random double between low and high
    static inline double randDouble(double low, double high) {
        return ((float)_int_distr(_eng) / numeric_limits<int>::max()) * (high - low) + low;
    }
    ///Get a random integer between low and high
    static inline int randInteger(int low, int high) {
        return low + (_int_distr(_eng) % (high - low + 1));
    }
};


#endif
