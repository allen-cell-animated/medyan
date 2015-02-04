
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

#ifndef M3SYM_ConjugateGradient_h
#define M3SYM_ConjugateGradient_h

#include "common.h"

#include "CGFletcherRievesMethod.h"
#include "CGPolakRibiereMethod.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

/// An implementation of [Minimzer](@ref Minimizer).
template <class CGType>
class ConjugateGradient : public Minimizer {
    
private:
    CGType _CGType;  ///< Implementation of a CG method
    double _GRADTOL; ///< Gradient tolerance used;
    
public:
    /// Constructor sets gradient tolerance parameter
    ConjugateGradient(double gradientTolerance) : _GRADTOL(gradientTolerance) {}
    
    void equlibrate(ForceFieldManager &FFM) {_CGType.minimize(FFM, _GRADTOL);}
};

#endif
