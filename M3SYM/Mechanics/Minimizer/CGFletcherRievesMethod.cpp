
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "CGFletcherRievesMethod.h"

#include "ForceFieldManager.h"
#include "Output.h"

void FletcherRieves::minimize(ForceFieldManager &FFM)
{
    
    int SpaceSize = 3 * BeadDB::instance()->size(); ///!!!!!! need to know
	double curEnergy = FFM.computeEnergy(0.0);
    cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.computeForces();

    //compute first gradient
	double gSquare = gradSquare();
    
	int numIter = 0;
	do {
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;

        //compute lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM);
        moveBeads(lambda);
        
        //compute new forces
        FFM.computeForcesAux();
        
        //compute new direction
		newGradSquare = gradAuxSquare();
		
        //calculate beta
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gSquare == 0) beta = 0;
            else beta = min(newGradSquare / gSquare, 1.0);
            
        }
        //shift gradient by beta
        if(gradDotProduct() < 0) shiftGradient(0.0);
        else shiftGradient(beta);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.computeEnergy(0.0);
        
		gSquare = newGradSquare;
        
	}
	while (gSquare > GRADTOL && (curEnergy - prevEnergy) < -ENERGYTOL);
    
	//cout << "Fletcher-Rieves Method: " << endl;
    //cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	//printForces();
	
}


