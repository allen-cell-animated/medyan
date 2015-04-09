
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
                                                       double MAXDIST)
{
    //system size
    int n = BeadDB::instance()->size();
    int ndof = 3 * n;
    if (ndof == 0) return;
    
    double curEnergy = FFM.computeEnergy(0.0);
    double prevEnergy;
    
    FFM.computeForces();
    
    //compute first gradient
    double curGrad = CGMethod::allFDotF();
    
    int numIter = 0;
    do {
        numIter++;
        double lambda, newGrad;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM, MAXDIST);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        
        //shift gradient
        shiftGradient(0.0);
        
        prevEnergy = curEnergy;
        curEnergy = FFM.computeEnergy(0.0);
        
        curGrad = newGrad;
    }
    while (curGrad / n > GRADTOL);
}

