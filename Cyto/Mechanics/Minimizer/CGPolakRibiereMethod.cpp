//
//  CGPolakRibiereMethod.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Output.h"

void PolakRibiere::minimize(ForceFieldManager &FFM){
    
    //cout<<"Forces before minimization:" <<endl;
	//PrintForces();
    
    int SpaceSize = 3 * BeadDB::instance()->size(); //// !!! change
	double curEnergy = FFM.computeEnergy(0.0);
    cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.computeForces();
    //printForces();
    
    //PrintForces();
	double gSquare = gradSquare();
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;
        
        lambda = backtrackingLineSearch(FFM);
        if(lambda < 0) lambda = 0.0;
        
        //cout<<"lambda= "<<lambda<<endl;
		//PrintForces();
        
        moveBeads(lambda);
        //PrintForces();
        
        FFM.computeForcesAux();
        //PrintForces();
        
		newGradSquare = gradAuxSquare();

		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gSquare == 0) beta = 0;
			else beta = max(0.0, (newGradSquare - gradDotProduct()/ gSquare));
        }
        //cout << "beta = " << beta <<endl;
		shiftGradient(beta);
        if(gradDotProduct() <= 0.0) shiftGradient(0);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.computeEnergy(0.0); 
        
        //PrintForces();
		gSquare = newGradSquare;
        //cout<<"GradSq before end=  "<< gSquare <<endl;
	}
	while (gSquare > GRADTOL && _energyChangeCounter <= ENERGYCHANGEITER);
    
	//cout << "Polak-Ribiere Method: " << endl;
    //cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    //printForces();
}
