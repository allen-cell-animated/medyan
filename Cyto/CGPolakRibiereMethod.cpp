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
		else
			beta = max(0.0, (newGradSquare - GradDotProduct()/ gradSquare));
		ShiftGradient(beta);
        
		prevVal = curVal;
		curVal = FFM.ComputeEnergy(0.0); /// Change maybe it to just compute energy and update energy or compute energyAux
        
		gradSquare = newGradSquare;
	}
	while (gradSquare > EPS);
    
    
	std::cout << "Polak-Ribiere Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    PrintForces();
}
