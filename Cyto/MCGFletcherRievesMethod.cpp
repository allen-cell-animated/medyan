//
//  MCGFletcherRievesMethod.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MCGFletcherRievesMethod.h"

using namespace std;
void FletcherRieves::Minimize(ForceFieldManager &FFM)
{
	
	const double EPS = 1e-10;
	
    int SpaceSize = 3 * BeadDB::Instance(getBeadDBKey())->size(); ///!!!!!! need to know
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

        lambda =0.01;
        //GoldenSection(pf);
        cout<<"lambda= "<<lambda<<endl;
		PrintForces();
        MoveBeads(lambda);
        PrintForces();
        
        FFM.ComputeForcesAux();
        PrintForces();
        
		newGradSquare = GradSquare(1);
		
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else
			beta = newGradSquare / gradSquare;
		ShiftGradient(beta);
        
		prevVal = curVal;
		curVal = FFM.ComputeEnergy(0.0);
        
		gradSquare = newGradSquare;
	}
	while (gradSquare > EPS);
    
    
	std::cout << "Fletcher-Rieves Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	PrintForces();
	
    
}


