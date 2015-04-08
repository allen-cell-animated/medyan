
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

void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL){

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
		double lambda, beta, allFADotFA, allFDotFA;
		vector<double> newGrad;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
        allFADotFA = CGMethod::allFADotFA();
        allFDotFA = CGMethod::allFDotFA();

        //choose beta
        //reset after ndof iterations
		if (numIter % ndof == 0)  beta = 0.0;
        
        //Polak-Ribieri update
        else beta = max(0.0, (allFADotFA - allFDotFA)/ curGrad);
    
        //shift gradient
        shiftGradient(0.0);
        
        //reset if direction not downhill
        if(CGMethod::allFDotFA() <= 0)
            shiftGradient(0.0);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.computeEnergy(0.0);
        
        curGrad = allFADotFA;
	}
	while (curGrad / n > GRADTOL);
}
