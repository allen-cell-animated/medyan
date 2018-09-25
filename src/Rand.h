
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
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

public:
    static mt19937 _eng;
    static long counter;
    static long Dcounter;
    static long Ncounter;
    ///Get a random double between low and high
    static inline double randDouble(double low, double high) {
        counter++;
        Dcounter++;
        return ((float)_int_distr(_eng) / numeric_limits<int>::max()) * (high - low) + low;
    }
    ///Get a random integer between low and high
    static inline int randInteger(int low, int high) {
        counter++;
        int y =_int_distr(_eng); 
        int x = low + (y % (high - low + 1));
        //std::cout<<"RandomInteger "<<x<<" "<<y<<" "<<high<<" "<<low<<" "<<counter<<" "<<Dcounter<<" "<<Ncounter<<endl;
        return x;
    }
};

#endif
