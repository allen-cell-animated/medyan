//
//  CGFletcherRievesMethod.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

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
    
    //PrintForces();
    
	double gSquare = gradSquare();
    //cout<<"GradSq=  "<<gradSquare<<endl;
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;

        lambda = backtrackingLineSearch(FFM);
        if(lambda < 0) return;
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
            else beta = newGradSquare / gSquare;
            
        }
		shiftGradient(beta);
        if(gradDotProduct() <= 0.0) shiftGradient(0);
        
		prevEnergy = curEnergy;
		curEnergy = FFM.computeEnergy(0.0);
        
		gSquare = newGradSquare;
        //cout<<"GradSq=  "<<gradSquare<<endl;
        
	}
	while (gSquare > GRADTOL && _energyChangeCounter <= ENERGYCHANGEITER);
    
	//cout << "Fletcher-Rieves Method: " << endl;
    //cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	//printForces();
	
}


