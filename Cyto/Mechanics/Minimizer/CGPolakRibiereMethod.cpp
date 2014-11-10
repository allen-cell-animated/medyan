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
    //Output o("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/beadoutput.txt");
    Output o("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/beadoutput.txt");
    o.printBasicSnapshot(0);
    
    int SpaceSize = 3 * BeadDB::instance(getBeadDBKey())->size(); //// !!! change
	double curEnergy = FFM.ComputeEnergy(0.0);
    //cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.ComputeForces();
    printForces();
    
    //PrintForces();
	double gSquare = gradSquare();
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;
        
        lambda = backtrackingLineSearch(FFM);
        if(lambda < 0) {
            cout<<"Lambda < 0" <<endl;
          break;
        }
        
        //cout<<"lambda= "<<lambda<<endl;
		//PrintForces();
        
        moveBeads(lambda);
        o.printBasicSnapshot(numIter);
        //PrintForces();
        
        FFM.ComputeForcesAux();
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
		curEnergy = FFM.ComputeEnergy(0.0); /// Change maybe it to just compute energy and update energy or compute energyAux
        
        //PrintForces();
		gSquare = newGradSquare;
        cout<<"GradSq before end=  "<< gSquare <<endl;
        cout << "Energy = " << curEnergy << endl;
        cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	}
	while (gSquare > GRADTOL && _energyChangeCounter <= ENERGYCHANGEITER);
    
	cout << "Polak-Ribiere Method: " << endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    printForces();
}
