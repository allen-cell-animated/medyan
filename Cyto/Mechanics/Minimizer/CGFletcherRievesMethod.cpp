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
	const double EPS = 1e-5;
	
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

        ///bracketing
        double ax = 0, bx = 0.0001, cx, fa, fb, fc;
        makeBracket(FFM, ax, bx, cx, fa, fb, fc);
        
        std::cout << "Bracket chosen: ax = " << ax << ", bx = " << bx << ", cx = "<< cx << std::endl;
        
        lambda = GoldenSection1(FFM, 1e-6);
        cout<<"lambda= "<<lambda<<endl;
        
		//PrintForces();
        MoveBeads(lambda);
        //PrintForces();
        //o.printBasicSnapshot(numIter);
        
        FFM.ComputeForcesAux();
        //PrintForces();
        
		newGradSquare = GradSquare(1);
		
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gradSquare == 0) beta = 0;
            else beta = newGradSquare / gradSquare;
            
        }
		ShiftGradient(beta);
        
		prevVal = curVal;
		curVal = FFM.ComputeEnergy(0.0);
        
		gradSquare = newGradSquare;
        cout<<"GradSq=  "<<gradSquare<<endl;
        
	}
	while (gradSquare > EPS);
    
    
	std::cout << "Fletcher-Rieves Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
	PrintForces();
	
    
}


