
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

void FletcherRieves::minimize(ForceFieldManager &FFM, double GRADTOL,
                                                      double ENERGYTOL,
                                                      double MAXDIST)
{
    //system size
    int N = Bead::numBeads();
    int NDOF = 3 * N;
    if (NDOF == 0) return;
    
    double curEnergy = FFM.computeEnergy(0.0);
    double prevEnergy;
    
    FFM.computeForces(); setBeads();
    
    //compute first gradient
    double curGrad = CGMethod::allFDotF();
    
    int numIter = 0;
    do {
        numIter++;
        double lambda, beta, newGrad;
        
        //find lambda by line search, move beads
        lambda = quadraticLineSearch(FFM, MAXDIST);
        moveBeads(lambda); setBeads();
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        
        //choose beta
        //reset after ndof iterations
        if (numIter % NDOF == 0)  beta = 0.0;
        //Fletcher-Rieves update
        else beta = newGrad / curGrad;
        
        //shift gradient
        shiftGradient(beta);
        
        //reset if search direction not downhill
        if(CGMethod::allFDotFA() <= 0)
            shiftGradient(0.0);
        
        prevEnergy = curEnergy;
        curEnergy = FFM.computeEnergy(0.0);
        
        curGrad = newGrad;
    }
    while (/* Iteration criterion */  numIter < 2 * NDOF &&
           /* Gradient tolerance  */  maxF() > GRADTOL &&
           /* Energy tolerance    */  curEnergy - prevEnergy <= -ENERGYTOL);
}
