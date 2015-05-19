
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


double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST) {
    
    double proj = allFDotFA();
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
        
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * proj;
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || idealEnergyChange >= -LSENERGYTOL)
            return 0.0;
    }
}

double CGMethod::quadraticLineSearch(ForceFieldManager& FFM, double MAXDIST) {
    
    double proj = allFDotFA();
    double f = maxF();
    
    //return zero if no forces
    if(f == 0.0) return 0.0;
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / f);
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //more vars
    double projOrig = proj;
    double projPrev = proj;
    double lambdaPrev = lambda;
    double energyPrevLambda = currentEnergy;
    
    double relErr, lambda0, delProj;
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        
        //compute new projection
        moveBeads(lambda);
        FFM.computeForcesAux();
        
        //move beads back
        resetBeads();
        
        proj = allFDotFA();
        delProj = proj - projPrev;
        
        if(fabs(proj) < EPS_QUAD || fabs(delProj) < EPS_QUAD)
            return 0.0;
            
        //check if ready for a quadratic projection
        relErr = fabs(1.0 - (0.5 * (lambda - lambdaPrev) * (proj + projPrev)+
                                           energyLambda) / energyPrevLambda);
        lambda0 = lambda - (lambda - lambdaPrev) * proj / delProj;
        
        //check if energy is decreasing and lambda within bounds
        if(relErr <= QUADRATICTOL &&
           energyLambda - currentEnergy <= 0 &&
           lambda0 > 0.0 && lambda0 < LAMBDAMAX)
            
            return lambda0;
    
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * projOrig;
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //save state
        projPrev = proj; lambdaPrev = lambda;
        energyPrevLambda = energyLambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || idealEnergyChange >= -LSENERGYTOL)
            return 0.0;
    }
}

