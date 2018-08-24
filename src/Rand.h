
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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

#include <limits>
#include <random>

/// A random number generator class.
class Rand {
    
private:
    static std::uniform_int_distribution<int> _int_distr;
    
public:
    static std::mt19937 eng; ///< The global random generator.
    static std::mt19937 engFixed; ///< Currently only used in tests.

    static long counter;
    static long Dcounter;
    static long Ncounter;

    ///Get a random double between low and high
    static inline double randDouble(double low, double high) {
        counter++;
        Dcounter++;
        return ((float)_int_distr(eng) / numeric_limits<int>::max()) * (high - low) + low;
    }
    ///Get a random integer between low and high
    static inline int randInteger(int low, int high) {
        counter++;
        int y =_int_distr(eng); 
        int x = low + (y % (high - low + 1));
        return x;
    }
};

#endif
