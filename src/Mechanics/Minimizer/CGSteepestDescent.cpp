
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

#include "CGSteepestDescent.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"

void SteepestDescent::minimize(ForceFieldManager &FFM, double GRADTOL,
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
    
    // TODO: update geometry here w/ derivative calculation
    FFM.computeForces();
    startMinimization();
    
    int numIter = 0;
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  maxF() > GRADTOL) {
        
        numIter++;
        double lambda;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);
        moveBeads(lambda); setBeads();
        
        // Update all geometry
        // TODO: update geometry here w/ derivative calculation
        
        // The unstretched local geometries have been calculated,
        // so forces could be calculated safely
        
        //compute new forces
        FFM.computeForcesAux();
        
        //shift gradient
        shiftGradient(0.0);
    }
    
    if (numIter >= N) {
        cout << endl;
        
        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->getParent()->printSelf();
        
        cout << "System energy..." << endl;
        FFM.computeEnergy(0.0, true);
        
        cout << endl;
    }
    
    //final force calculation
    FFM.computeForces();
    FFM.computeLoadForces();
}

