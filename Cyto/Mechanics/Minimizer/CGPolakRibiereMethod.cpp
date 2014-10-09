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
    const double EPS = 1e-4;
    //Output o("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/beadoutput.txt");
    //o.printBasicSnapshot(0);
    
    int SpaceSize = 3 * BeadDB::Instance(getBeadDBKey())->size(); //// !!! change
	double curEnergy = FFM.ComputeEnergy(0.0);
    //cout<<"Energy = "<< curEnergy <<endl;
	double prevEnergy = curEnergy;
	FFM.ComputeForces();
    
    //PrintForces();
    
    _energyChangeCounter = 0;
	double gradSquare = GradSquare();
    //cout<<"GradSq=  "<<gradSquare<<endl;
    
	int numIter = 0;
	do
	{
		numIter++;
		double lambda, beta, newGradSquare;
		vector<double> newGrad;
        
        ///initial bracketing
        //double ax = 0, bx =0.0001 , cx=100, fa, fb, fc;
        //makeBracket(FFM, ax, bx, cx, fa, fb, fc);
        
        //std::cout << "Bracket chosen: ax = " << ax << ", bx = " << bx << ", cx = "<< cx << std::endl;
        
        lambda = BacktrackingLineSearch(FFM);
        if(lambda < 0) {
            std::cout << "Line search stopping." <<std::endl;
            break;
        }
//        if(lambda == 0) {
//            std::cout << "Lambda is zero." << std::endl;
//            //break;
//        }
        
        //cout<<"lambda= "<<lambda<<endl;
		//PrintForces();
        
        //cout<<"GradSq before move beads=  "<<gradSquare<<endl;
        
        MoveBeads(lambda);
//        if (lambda > _lambdaMin) EnergyBacktracking(FFM, lambda, curEnergy);
        
        //o.printBasicSnapshot(numIter);
        //PrintForces();
        
        FFM.ComputeForcesAux();
        //PrintForces();
        
		newGradSquare = GradSquare(1);
		
        //cout << "Grad dot prod "<< GradDotProduct() << endl;
        //cout << "New grad square = " << newGradSquare << endl;
		if (numIter % (5 * SpaceSize) == 0) beta = 0;
		else {
            if(gradSquare == 0) beta = 0;
			else beta = max(0.0, (newGradSquare - GradDotProduct()/ gradSquare));
        }
        //cout << "beta = " << beta <<endl;
		ShiftGradient(beta);
        if(GradDotProduct() <= 0.0) ShiftGradient(0);
        
		prevEnergy = curEnergy;
        //cout << "Calling last compute energy in minimizer" << endl;
		curEnergy = FFM.ComputeEnergy(0.0); /// Change maybe it to just compute energy and update energy or compute energyAux
        
        //PrintForces();
		gradSquare = newGradSquare;
        //cout<<"GradSq before end=  "<<gradSquare<<endl;
        //cout << "Energy = " << curEnergy << endl;
        //cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
        
	}
	while (gradSquare > EPS && _energyChangeCounter <= _maxEnergyChangeIter);
    
	std::cout << "Polak-Ribiere Method: " << std::endl;
    cout<<"numIter= " <<numIter<<"  Spacesize = "<<SpaceSize <<endl;
    PrintForces();
}
