
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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Output.h"

void PolakRibiere::minimize(ForceFieldManager &FFM){
    
    int SpaceSize = 3 * BeadDB::instance()->size(); //// !!! change
	double curEnergy = FFM.computeEnergy(0.0);
    cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.computeForces();

    //compute first gradient
	double gSquare = gradSquare();
    
	int numIter = 0;
	do {
		numIter++;
		double lambda, beta, newGradSquare, conjSquare;
		vector<double> newGrad;
        
        //find lambda by line search, move beads
        lambda = backtrackingLineSearch(FFM);
        //cout << "Lamba" << lambda << endl;
        moveBeads(lambda);

        //compute new forces
        FFM.computeForcesAux();
        
        //compute direction
		newGradSquare = gradAuxSquare();
        conjSquare = gradDotProduct();

        //choose beta, safeguard for blowups
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gSquare == 0) beta = 0;
            else beta = min(max(0.0, (newGradSquare - conjSquare)/ gSquare), 1.0);
        }
        //cout << "Beta = " << beta << endl;
        shiftGradient(beta);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.computeEnergy(0.0); 
		gSquare = newGradSquare;
        cout << "Current energy = " << curEnergy << endl;
        //cout << "Previous energy = " << prevEnergy << endl;
        //cout << "GradSquare = " << gSquare << endl;
	}
	while (gSquare > GRADTOL && (curEnergy - prevEnergy) < -ENERGYTOL);

    cout<<"Energy = "<< curEnergy <<endl;
//    
//	cout << "Polak-Ribiere Method: " << endl;
//  cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
//  printForces();
}
