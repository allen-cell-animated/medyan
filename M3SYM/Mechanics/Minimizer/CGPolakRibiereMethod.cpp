
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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Output.h"

void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                                                    double MAXDIST,
                                                    double LAMBDAMAX){
    
    //system size
    int N = Bead::numBeads();
    
	FFM.computeForces();
    startMinimization();

    //compute first gradient
    double curGrad = CGMethod::allFDotF();
    
	int numIter = 0;
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  maxF() > GRADTOL) {

		numIter++;
		double lambda, beta, newGrad, prevGrad;
        
        //find lambda by line search, move beads
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);
        
        moveBeads(lambda); setBeads();
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();
        
        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        
        //update prev forces
        FFM.computeForcesAuxP();
        
        //shift gradient
        shiftGradient(beta);
        
        //direction reset if not downhill, or no progress made
        if(CGMethod::allFDotFA() <= 0 || areSame(curGrad, newGrad)) {
            
            shiftGradient(0.0);
            _safeMode = true;
        }
        
        curGrad = newGrad;
    }
    
    if (numIter >= N) {
        cout << "WARNING: Did not minimize in N (= number of beads) steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;
        
        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->printInfo();
    }
}
