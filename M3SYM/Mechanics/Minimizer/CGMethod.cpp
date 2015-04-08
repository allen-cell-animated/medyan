
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
	for(auto b: *BeadDB::instance())
        g += b->FDotF();
    return g;
}

double CGMethod::allFADotFA()
{
    double g = 0;
	for(auto b: *BeadDB::instance())
        g += b->FADotFA();
    
    return g;
}

double CGMethod::allFDotFA()
{
    double g = 0;
	for(auto b: *BeadDB::instance())
        g += b->FDotFA();
    
    return g;
}


void CGMethod::moveBeads(double d)
{
	for(auto b: *BeadDB::instance()) {
        
        b->coordinate[0] = b->coordinate[0] + d* b->force[0];
        b->coordinate[1] = b->coordinate[1] + d* b->force[1];
        b->coordinate[2] = b->coordinate[2] + d* b->force[2];
        
        //reset coord Aux
        b->coordinateAux = b->coordinate;
	}
}

void CGMethod::moveBeadsAux(double d) {
    
    for(auto b: *BeadDB::instance()) {
        
        b->coordinateAux[0] = b->coordinate[0] + d* b->force[0];
        b->coordinateAux[1] = b->coordinate[1] + d* b->force[1];
        b->coordinateAux[2] = b->coordinate[2] + d* b->force[2];
	}
}


void CGMethod::shiftGradient(double d)
{
	for(auto b: *BeadDB::instance()) {
        
        b->force[0] = b->forceAux[0] + d* b->force[0];
        b->force[1] = b->forceAux[1] + d* b->force[1];
        b->force[2] = b->forceAux[2] + d* b->force[2];
	}
}

void CGMethod::printForces()
{
	cout << "Print Forces" << endl;
    for(auto b: *BeadDB::instance()) {
        
		for (int i = 0; i<3; i++)
            cout << b->coordinate[i] << "  "<<
                    b->force[i] <<"  "<<b->forceAux[i]<<endl;
	}
    cout << "End of Print Forces" << endl;
}


double CGMethod::goldenSection(ForceFieldManager& FFM)
{
	double a = 0;
	double b = 200;
    double phi = 0.5 * (1 + sqrt(5) );
    double inv_phi = 1/phi;
	double x1 = b - inv_phi * (b - a);
	double x2 = a + inv_phi * (b - a);
    
	while (fabs(b - a) > LSENERGYTOL)
	{
		if (FFM.computeEnergy(x1) >= FFM.computeEnergy(x2) ){
            a = x1;
            x1 = x2;
            x2 =a + inv_phi * (b - a);
        }
        else {
            b = x2;
            x2 = x1;
            x1 = b - inv_phi * (b - a);
        }
    }
    double returnLambda = (a + b)/2.0;
    
    return returnLambda;
}


double CGMethod::binarySearch(ForceFieldManager& FFM)
{
    double a = 0, b = 100;
    while (fabs(b - a) > LSENERGYTOL){
        
        double half_x1 = ((a + b)/2 - LSENERGYTOL/4);
        double half_x2 = ((a + b)/2 + LSENERGYTOL/4);
        if (FFM.computeEnergy(half_x1) <= FFM.computeEnergy(half_x2)) b = half_x2;
        else a = half_x1;
    }
     return (a + b) / 2;
}


double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM) {
    
    double allFDotFA = 0.0; double maxF = 0.0;
    
    for(auto b: *BeadDB::instance()) {
        
        //calc f dot fa
        allFDotFA += b->FDotFA();
        
        //calc max force
        for(int i = 0 ; i < 3; i++)
            maxF = max(maxF, fabs(b->force[i]));
    }
    
    //return zero if no forces
    if(maxF == 0.0) return 0.0;
 
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / maxF);
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA;
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange)
            return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || idealEnergyChange >= -LSENERGYTOL)
            return 0.0;
    }
}

double CGMethod::quadraticLineSearch(ForceFieldManager& FFM) {
    
    double allFDotFA = 0.0; double maxF = 0.0;
    
    for(auto b: *BeadDB::instance()) {
        
        //calc f dot fa
        allFDotFA += b->FDotFA();
        
        //calc max force
        for(int i = 0 ; i < 3; i++)
            maxF = max(maxF, fabs(b->force[i]));
    }
    
    //return zero if no forces
    if(maxF == 0.0) return 0.0;
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / maxF);

    //tracking conjugate direction
    double allFDotFAPrev = allFDotFA;
    double allFDotFAOrig = allFDotFA;
    
    double deltaFDotFA, relErr, lambda0;
    
    double energyOrig = FFM.computeEnergy(0.0);
    
    double energyLambda;
    double energyPrevLambda = energyOrig;
    double idealEnergyChange, energyChange;
    
    double lambdaPrev = 0.0;
    
    //backtracking loop
    while(true) {
        //move beads
        energyLambda = FFM.computeEnergy(lambda);
        moveBeadsAux(lambda);
        
        //calculate new gradients
        FFM.computeForcesAux();
        allFDotFA = CGMethod::allFDotFA();
        
        //compute gradient change
        deltaFDotFA = allFDotFA - allFDotFAPrev;
        
        //if no gradient change, return 0.0
        if(fabs(allFDotFA) < EPS_QUAD || fabs(deltaFDotFA) < EPS_QUAD)
            return 0.0;
        
        //Check if ready for a quadratic projection
        if(deltaFDotFA <= 0) {
        
            relErr = fabs(1.0 - (0.5 * (lambda - lambdaPrev) *
                                 (allFDotFA + allFDotFAPrev) +
                                  energyLambda) / energyPrevLambda);
            
            lambda0 = lambda - (lambda - lambdaPrev) *
                      allFDotFA / deltaFDotFA;
            
            if(relErr <= QUADRATICTOL &&
               0.0 < lambda0 && lambda0 < LAMBDAMAX)
                return lambda0;
        }
        
        //calculate ideal change
        idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFAOrig;
        energyChange = energyLambda - energyOrig;
        
        //return if ok
        if(energyChange <= idealEnergyChange)
            return lambda;
            
        //save previous state
        allFDotFAPrev = allFDotFA;
        energyPrevLambda = energyLambda;
        
        //reduce lambda
        lambdaPrev = lambda;
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || idealEnergyChange >= -LSENERGYTOL)
            return 0.0;
    }
}

