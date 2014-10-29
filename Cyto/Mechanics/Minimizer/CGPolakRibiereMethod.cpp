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

using namespace std;
void PolakRibiere::Minimize(ForceFieldManager &FFM){
    
    //cout<<"Forces before minimization:" <<endl;
	//PrintForces();
    Output o("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/beadoutput.txt");
    o.printBasicSnapshot(0);
    
    int SpaceSize = 3 * BeadDB::Instance(getBeadDBKey())->size(); //// !!! change
	double curEnergy = FFM.ComputeEnergy(0.0);
    //cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.ComputeForces();
    
    //PrintForces();
	double gradSquare = GradSquare();
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;
        
        lambda = QuadraticLineSearch(FFM);
        if(lambda < 0) {
            cout<<"Lambda < 0" <<endl;
          break;
        }
        
        cout<<"lambda= "<<lambda<<endl;
		//PrintForces();
        
        MoveBeads(lambda);
        o.printBasicSnapshot(numIter);
        //PrintForces();
        
        FFM.ComputeForcesAux();
        //PrintForces();
        
		newGradSquare = GradAuxSquare();

		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gradSquare == 0) beta = 0;
			else beta = max(0.0, (newGradSquare - GradDotProduct()/ gradSquare));
        }
        //cout << "beta = " << beta <<endl;
		ShiftGradient(beta);
        if(GradDotProduct() <= 0.0) ShiftGradient(0);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.ComputeEnergy(0.0); /// Change maybe it to just compute energy and update energy or compute energyAux
        
        //PrintForces();
		gradSquare = newGradSquare;
        cout<<"GradSq before end=  "<<gradSquare<<endl;
        cout << "Energy = " << curEnergy << endl;
        cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
        
	}
	while (gradSquare > GRADTOL && _energyChangeCounter <= ENERGYCHANGEITER);
    
	std::cout << "Polak-Ribiere Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    PrintForces();
}
