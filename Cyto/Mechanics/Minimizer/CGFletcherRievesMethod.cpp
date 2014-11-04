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

using namespace std;
void FletcherRieves::Minimize(ForceFieldManager &FFM)
{
    //Output o("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/beadoutput.txt");
    //o.printBasicSnapshot(0);
    
    int SpaceSize = 3 * BeadDB::Instance(getBeadDBKey())->size(); ///!!!!!! need to know
	double curVal = FFM.ComputeEnergy(0.0);
    //cout<<"Energy = "<< curVal <<endl;
	double prevVal = curVal;
	FFM.ComputeForces();
    
    //PrintForces();
    
	double gradSquare = GradSquare();
    //cout<<"GradSq=  "<<gradSquare<<endl;
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;

        lambda = BacktrackingLineSearch(FFM);
        if(lambda < 0) break;
        //cout<<"lambda= "<<lambda<<endl;
        
		//PrintForces();
        MoveBeads(lambda);
        //PrintForces();
        //o.printBasicSnapshot(numIter);
        
        FFM.ComputeForcesAux();
        //PrintForces();
        
		newGradSquare = GradAuxSquare();
		
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gradSquare == 0) beta = 0;
            else beta = newGradSquare / gradSquare;
            
        }
		ShiftGradient(beta);
        if(GradDotProduct() <= 0.0) ShiftGradient(0);
        
		prevVal = curVal;
		curVal = FFM.ComputeEnergy(0.0);
        
		gradSquare = newGradSquare;
        //cout<<"GradSq=  "<<gradSquare<<endl;
        
	}
	while (gradSquare > GRADTOL && _energyChangeCounter <= ENERGYCHANGEITER);
    
	std::cout << "Fletcher-Rieves Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	PrintForces();
	
}


