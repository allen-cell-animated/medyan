
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

#ifndef MEDYAN_ConjugateGradient_h
#define MEDYAN_ConjugateGradient_h

#include "common.h"

#include "CGFletcherRievesMethod.h"
#include "CGPolakRibiereMethod.h"
#include "CGSteepestDescent.h"
#include "Minimizer.h"

#include "ForceFieldManager.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

/// An implementation of [Minimzer](@ref Minimizer).
template <class CGType> class ConjugateGradient : public Minimizer {
    
private:
    CGType _CGType;  ///< Implementation of a CG method
    
    floatingpoint _GRADTOL;   ///< Gradient tolerance used
    floatingpoint _MAXDIST;   ///< Max distance used to move
    floatingpoint _LAMBDAMAX; ///< Maximum lambda that can be returned
    floatingpoint _LAMBDARUNNINGAVERAGEPROBABILITY = 0.0; //running average probability.
    
public:
    /// Constructor sets gradient tolerance parameter
    ConjugateGradient(floatingpoint gradientTolerance,
                      floatingpoint maxDistance,
                      floatingpoint lambdaMax,
                      floatingpoint lambdarunningaverageprobability)
    
        : _GRADTOL(gradientTolerance),
          _MAXDIST(maxDistance),
          _LAMBDAMAX(lambdaMax),
          _LAMBDARUNNINGAVERAGEPROBABILITY(lambdarunningaverageprobability) {}
    
    
    /// This function will minimize the system until the following criterion are met:
    /// 1) Largest force in the network < GRADTOL
    /// 3) Number of iterations exceeds 5N, unless in initial minimization
    void equlibrate(ForceFieldManager &FFM, bool steplimit) {
        _CGType.minimize(FFM, _GRADTOL, _MAXDIST, _LAMBDAMAX,
                _LAMBDARUNNINGAVERAGEPROBABILITY, steplimit);
    }





    tuple<floatingpoint, vector<floatingpoint>, vector<string>> getEnergy(ForceFieldManager &FFM, floatingpoint d){
      
        //double* coord = _CGType.getCoords();
        floatingpoint* coord = Bead::getDbData().coords.data();
        
        FFM.vectorizeAllForceFields();

        tuple<floatingpoint, vector<floatingpoint>, vector<string>> HRMDvec = FFM.computeEnergyHRMD(coord);
        
        // delete [] coord;
        
        FFM.cleanupAllForceFields();
        
        
        
        return HRMDvec;
        
        
        
    }

    
    
};

#endif
