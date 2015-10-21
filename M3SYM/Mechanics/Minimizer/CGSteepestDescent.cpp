
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

#include "CGSteepestDescent.h"

#include "ForceFieldManager.h"
#include "Output.h"

void SteepestDescent::minimize(ForceFieldManager &FFM, double GRADTOL,
                                                       double MAXDIST,
                                                       double LAMBDAMAX)
{
    //system size
    int N = Bead::numBeads();
    
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
        
        //compute new forces
        FFM.computeForcesAux();
        
        //shift gradient
        shiftGradient(0.0);
    }
    
    if (numIter >= N) {
        cout << "WARNING: Did not minimize in N (= number of beads) steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        maxBead()->printSelf();
        
        cout << "System energy..." << endl;
        FFM.computeEnergy(0.0, true);
        
        cout << endl;
    }
}

