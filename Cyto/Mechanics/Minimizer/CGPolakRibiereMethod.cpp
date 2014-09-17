//
//  CGPolakRibiereMethod.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CGPolakRibiereMethod.h"

using namespace std;
void PolakRibiere::Minimize(ForceFieldManager &FFM){
    
    cout<<"Forces before minimization:" <<endl;
	PrintForces();
    const double EPS = 1e-10;
	
    int SpaceSize = 3 * BeadDB::Instance(getBeadDBKey())->size(); //// !!! change
	double curVal = FFM.ComputeEnergy(0.0);
    cout<<"Energy = "<< curVal <<endl;
	double prevVal = curVal;
	FFM.ComputeForces();
    
    PrintForces();
    
	double gradSquare = GradSquare();
    
    cout<<"GradSq=  "<<gradSquare<<endl;
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;
        
        lambda = GoldenSection(FFM);
        cout<<"lambda= "<<lambda<<endl;
		PrintForces();
        
        MoveBeads(lambda);
        
        PrintForces();
        
        FFM.ComputeForcesAux();
        PrintForces();
        
		newGradSquare = GradSquare(1);
		
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gradSquare == 0) beta = 0;
			else beta = max(0.0, (newGradSquare - GradDotProduct()/ gradSquare));
        }
		ShiftGradient(beta);
        
		prevVal = curVal;
        cout << "Calling last compute energy in minimizer" << endl;
		curVal = FFM.ComputeEnergy(0.0); /// Change maybe it to just compute energy and update energy or compute energyAux
        
        PrintForces();
		gradSquare = newGradSquare;
        cout<<"GradSq before end=  "<<gradSquare<<endl;

        
	}
	while (gradSquare > EPS);
    
    
	std::cout << "Polak-Ribiere Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    PrintForces();
}
