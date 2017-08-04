
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

#include "CGFletcherRievesMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"

void FletcherRieves::minimize(ForceFieldManager &FFM, double GRADTOL,
                              double MAXDIST, double LAMBDAMAX, bool steplimit)
{
    //number of steps
    int N;
    if(steplimit) {
        int beadMaxStep = 5 * Bead::numBeads();
        N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
    }
    else {
        N = numeric_limits<int>::max();
    }
    
    startMinimization();
    FFM.vectorizeAllForceFields();
    
    FFM.computeForces(coord, force);
    FFM.copyForces(forceAux, force);
    
    //compute first gradient
    double curGrad = CGMethod::allFDotF();
    
    int numIter = 0;
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  maxF() > GRADTOL) {
        numIter++;
        double lambda, beta, newGrad;
        
        //find lambda by line search, move beads
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForces(coord, forceAux);
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        
        //Fletcher-Rieves update
        beta = newGrad / curGrad;
        
        //shift gradient
        shiftGradient(beta);
        
        //direction reset if not downhill or no progress made
        if(CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)) {
            shiftGradient(0.0);
            _safeMode = true;
        }
        
        curGrad = newGrad;
    }
    
    if (numIter >= N) {
        cout << endl;
        
        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->getParent()->printSelf();
        
        cout << "System energy..." << endl;
        FFM.computeEnergy(coord, force, 0.0, true);
        
        cout << endl;
    }
    
    //final force calculation
    FFM.computeForces(coord, force);
    FFM.copyForces(forceAux, force);
    FFM.computeLoadForces();
    endMinimization();
    
    FFM.cleanupAllForceFields();
}
