
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

#include <cmath>
#include <limits>
#include <random>
#include <type_traits>

#include "common.h"

/// A random number generator class.
class Rand {
    
private:
    static uniform_int_distribution<int> _int_distr;

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
};

// Safe exponential distribution
//
// Generates exponential distribution safely:
//   - Precondition: λ is not negative (not checked)
//   - if λ is zero: return infinity
//   - if λ is positive: set the λ for the specified exponential distribution,
//     and generate a random number which is not infinity.
//
// Inputs:
//   - d:      the exponential distribution object
//   - lambda: the exponential distribution parameter
//   - g:      the random number generator
//
// Output:
//   - The generated random number with type ExpDist::result_type
//
// Note:
//   - When λ = 0, the input distribution and generator are not touched.
//   - When λ > 0, the generator can be used 1 or more times.
//
// Background:
//   - std::exponential_distribution requires that λ > 0, though most
//     implementations would return ∞ when λ = 0.
//   - Some implementations of std::exponential_distribution returns infinity
//     when λ is positive. See
//     http://open-std.org/JTC1/SC22/WG21/docs/lwg-active.html#2524
//-----------------------------------------------------------------------------
template< typename ExpDist, typename FloatLambda, typename Generator >
inline auto safeExpDist(
    ExpDist&    d,
    FloatLambda lambda,
    Generator&  g
) {
    using ResType = typename ExpDist::result_type;

    ResType res;
    if(lambda == 0) {
        res = std::numeric_limits< ResType >::infinity();
    } else {
        // Set d's parameter
        d.param(typename ExpDist::param_type( lambda ));

        // Generate non-inf random number
        do res = d(g); while(std::isinf(res));
    }

    return res;
}

#endif
