
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

#include "CGFletcherRievesMethod.h"

#include "ForceFieldManager.h"
#include "Output.h"

void FletcherRieves::minimize(ForceFieldManager &FFM, double GRADTOL)
{
    //system size
    int n = BeadDB::instance()->size();
    int ndof = 3 * n;
    if (ndof == 0) return;
    
    double curEnergy = FFM.computeEnergy(0.0);
    double prevEnergy;
    
    FFM.computeForces();
    
    //compute first gradient
    double curGradient = CGMethod::allFDotF();
    
    int numIter = 0;
    do {
        numIter++;
        double lambda, beta, newGradient;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        newGradient = CGMethod::allFADotFA();
        
        //choose beta
        //reset after ndof iterations
        if (numIter % ndof == 0)  beta = 0.0;
        
        //Fletcher-Rieves update
        else beta = newGradient / curGradient;
        
        //shift gradient
        shiftGradient(beta);
        
        prevEnergy = curEnergy;
        curEnergy = FFM.computeEnergy(0.0);
        
        curGradient = newGradient;
    }
    while (curGradient / n > GRADTOL);
}
