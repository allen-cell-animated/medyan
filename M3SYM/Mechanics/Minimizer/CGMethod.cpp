
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "CGMethod.h"

#include "ForceFieldManager.h"
#include "Bead.h"


double CGMethod::allFDotF()
{
    double g = 0;
    for(auto b: Bead::getBeads())
        g += b->FDotF();
    return g;
}

double CGMethod::allFADotFA()
{
    double g = 0;
	for(auto b: Bead::getBeads())
        g += b->FADotFA();
    
    return g;
}

double CGMethod::allFADotFAP()
{
    double g = 0;
    for(auto b: Bead::getBeads())
        g += b->FADotFAP();
    
    return g;
}

double CGMethod::allFDotFA()
{
    double g = 0;
	for(auto b: Bead::getBeads())
        g += b->FDotFA();
    
    return g;
}

double CGMethod::maxF() {
    
    double maxF = 0;
    
    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));

    return maxF;
}

void CGMethod::startMinimization() {

    for(auto b: Bead::getBeads())
        b->coordinateB = b->coordinateP = b->coordinate;
}


void CGMethod::moveBeads(double d)
{
	for(auto b: Bead::getBeads()) {
        
        b->coordinate[0] = b->coordinate[0] + d * b->force[0];
        b->coordinate[1] = b->coordinate[1] + d * b->force[1];
        b->coordinate[2] = b->coordinate[2] + d * b->force[2];
	}
}

void CGMethod::resetBeads() {
    
    for(auto b: Bead::getBeads())
        b->coordinate = b->coordinateP;
}

void CGMethod::setBeads() {
    
    for(auto b: Bead::getBeads())
        b->coordinateP = b->coordinate;
}


void CGMethod::shiftGradient(double d)
{
	for(auto b: Bead::getBeads()) {
        
        b->force[0] = b->forceAux[0] + d * b->force[0];
        b->force[1] = b->forceAux[1] + d * b->force[1];
        b->force[2] = b->forceAux[2] + d * b->force[2];
	}
}

void CGMethod::printForces()
{
	cout << "Print Forces" << endl;
    for(auto b: Bead::getBeads()) {
        
		for (int i = 0; i<3; i++)
            cout << b->coordinate[i] << "  "<<
                    b->force[i] <<"  "<<b->forceAux[i]<<endl;
	}
    cout << "End of Print Forces" << endl;
}

double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                double LAMBDAMAX) {
    
    double f = maxF();
    
    //return zero if no forces
    if(f == 0.0) return 0.0;
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / f);
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= 0.0) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || lambda <= LAMBDATOL) {
            return 0.0;
        }
    }
}
