
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
                                                      double MAXDIST)
{
    //system size
    int n = BeadDB::instance()->size();
    int ndof = 3 * n;
    if (ndof == 0) return;
    
    double curEnergy = FFM.computeEnergy(0.0);
    double prevEnergy;
    
    FFM.computeForces();
    
    cout << "Initial energy = " << curEnergy << endl;
    
    //compute first gradient
    double curGrad = CGMethod::allFDotF();
    
    int numIter = 0;
    do {
        numIter++;
        double lambda, beta, newGrad, prevGrad;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM, MAXDIST);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();
        
        //choose beta
        //reset after ndof iterations
        if (numIter % ndof == 0)  beta = 0.0;
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
    while (curGrad / n > GRADTOL);
}
