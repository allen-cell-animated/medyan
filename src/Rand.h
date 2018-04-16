
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

#ifndef MEDYAN_Rand_h
#define MEDYAN_Rand_h

#include <limits>
#include <random>

/// A random number generator class.
class Rand {
    
private:
    static std::uniform_int_distribution<int> _int_distr;
    
public:
    static std::mt19937 eng;
    static std::mt19937 engFixed;

    ///Get a random double between low and high
    static inline double randDouble(double low, double high) {
        return ((float)_int_distr(eng) / std::numeric_limits<int>::max()) * (high - low) + low;
    }
    ///Get a random integer between low and high
    static inline int randInteger(int low, int high) {
        return low + (_int_distr(eng) % (high - low + 1));
    }
};


#endif
