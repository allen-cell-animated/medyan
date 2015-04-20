
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
#include "CGSteepestDescent.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

/// An implementation of [Minimzer](@ref Minimizer).
template <class CGType>
class ConjugateGradient : public Minimizer {
    
private:
    CGType _CGType;  ///< Implementation of a CG method
    double _ENERGYTOL; ///< Energy tolerance used
    double _GRADTOL;   ///< Gradient tolerance used
    double _MAXDIST;   ///< Max distance used to move
    
public:
    /// Constructor sets gradient tolerance parameter
    ConjugateGradient(double gradientTolerance,
                      double energyTolerance,
                      double maxDistance)
    
        : _GRADTOL(gradientTolerance),
          _ENERGYTOL(energyTolerance),
          _MAXDIST(maxDistance) {}
    
    
    /// This function will minimize the system until the following criterion are met:
    /// 1) Sum of all forces in network are < GRADTOL
    /// 2) Difference in energy between steps < -ENERGYTOL
    /// 3) Number of iterations exceeds 2 * NDOF
    void equlibrate(ForceFieldManager &FFM) {
        _CGType.minimize(FFM, _GRADTOL, _ENERGYTOL, _MAXDIST);
    }
};

#endif
