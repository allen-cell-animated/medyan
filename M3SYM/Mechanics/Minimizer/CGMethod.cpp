
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "CGMethod.h"

#include "ForceFieldManager.h"
#include "Bead.h"


double CGMethod::gradSquare()
{
    double g = 0;
	for(auto it: *BeadDB::instance())
        g += (*it).calcForceSquare();
    return g;
}

double CGMethod::gradAuxSquare()
{
    double g = 0;
	for(auto it: *BeadDB::instance())
        g += (*it).calcForceAuxSquare();
    
    return g;
}

double CGMethod::gradDotProduct()
{
    double g = 0;
	for(auto it: *BeadDB::instance())
        g += (*it).calcDotForceProduct();
    
    return g;
}


void CGMethod::moveBeads(double d)
{
	for(auto it: *BeadDB::instance()) {
        
        (*it).coordinate[0] = (*it).coordinate[0] + d* (*it).force[0];
        (*it).coordinate[1] = (*it).coordinate[1] + d* (*it).force[1];
        (*it).coordinate[2] = (*it).coordinate[2] + d* (*it).force[2];
        
        //reset coord Aux
        (*it).coordinateAux = (*it).coordinate;
	}
}

void CGMethod::moveBeadsAux(double d) {
    
    for(auto it: *BeadDB::instance()) {
        
        (*it).coordinateAux[0] = (*it).coordinateAux[0] + d* (*it).force[0];
        (*it).coordinateAux[1] = (*it).coordinateAux[1] + d* (*it).force[1];
        (*it).coordinateAux[2] = (*it).coordinateAux[2] + d* (*it).force[2];
	}
}


void CGMethod::shiftGradient(double d)
{
	for(auto it: *BeadDB::instance()) {
        
        (*it).force[0] = (*it).forceAux[0] + d* (*it).force[0];
        (*it).force[1] = (*it).forceAux[1] + d* (*it).force[1];
        (*it).force[2] = (*it).forceAux[2] + d* (*it).force[2];
	}
}

void CGMethod::printForces()
{
	cout << "Print Forces" << endl;
    for(auto it: *BeadDB::instance()) {
        
		for (int i = 0; i<3; i++)
            cout << (*it).coordinate[i] << "  "<< (*it).force[i]<<"  "
                 <<(*it).forceAux[i]<<endl;
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
    
    //check if return value is in bounds of lambda min and max
    if (returnLambda > LAMBDAMAX) return LAMBDAMAX;
    else if(returnLambda < LAMBDAMIN) return LAMBDAMIN;
    else return returnLambda;
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
    
    double conjDirection = 0.0;
    double maxForce = 0.0;
    
    for(auto it: *BeadDB::instance()) {
        conjDirection += it->calcDotForceProduct();
        for(int i=0 ; i < 3; i++) maxForce = max(maxForce, fabs(it->force[i]));
    }
    //return zero if no forces
    if(maxForce == 0.0) return 0.0;
 
    //calculate first lambda. cannot be greater than
    //lambda max, less than lambdamin
    double lambda = min(LAMBDAMAX, MAXDIST / maxForce);
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        //cout << "EnergyLambda = " << energyLambda << endl;
        
        //calculate ideal change
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * conjDirection;
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        if(lambda <= LAMBDAMIN || idealEnergyChange >= -LSENERGYTOL) {
            
            if(energyChange < 0) {
                //see if we can return LAMBDAMIN
                if(FFM.computeEnergy(LAMBDAMIN) <= currentEnergy) return LAMBDAMIN;
                else return lambda;
            }
            //can't, just return 0
            else return 0.0;
        }
    }
}

double CGMethod::quadraticLineSearch(ForceFieldManager& FFM) {
    
    double conjDirection = 0.0;
    double maxForce = 0.0;
    
    for(auto it: *BeadDB::instance()) {
        conjDirection += it->calcDotForceProduct();
        for(int i=0 ; i < 3; i++) maxForce = max(maxForce, fabs(it->force[i]));
    }
    //return zero if no forces
    if(maxForce == 0.0) return 0.0;
    
    //set up backtracking
    double lambda = min(LAMBDAMAX, MAXDIST / maxForce);
    double lambdaPrev = 0.0;
    
    double conjDirectionPrev = conjDirection;
    double deltaconjDirection, relErr, lambda0;
    
    double energyCurrent = FFM.computeEnergy(0.0);
    double energyPrev = energyCurrent;
    double idealEnergyChange, energyChange;
    
    //backtracking loop
    while(true) {
        
        //move beads
        energyCurrent = FFM.computeEnergy(lambda);
        moveBeadsAux(lambda);
        
        //calculate new gradients
        FFM.computeForcesAux();
        conjDirection = gradDotProduct();
        
        //compute gradient change
        deltaconjDirection = conjDirection - conjDirectionPrev;
        
        //if no gradient change, return 0.0
        if(fabs(conjDirection) < EPS_QUAD || fabs(deltaconjDirection) < EPS_QUAD)
            return 0.0;
        
        //Check if ready for a quadratic projection
        relErr = fabs(1.0 - (0.5 * (lambda - lambdaPrev) *
                     (conjDirection + conjDirectionPrev) + energyCurrent) / energyPrev);
        lambda0 = lambda - (lambda - lambdaPrev) * conjDirection / deltaconjDirection;
        
        if(relErr <= QUADRATICTOL && lambda0 > 0.0 && lambda0 < LAMBDAMAX)
            return lambda0;
        
        //calculate ideal change
        idealEnergyChange = -BACKTRACKSLOPE * lambda * conjDirection;
        energyChange = energyCurrent - energyPrev;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //save previous state
        conjDirectionPrev = conjDirection;
        energyPrev = energyCurrent;
        lambdaPrev = lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        if(lambda <= LAMBDAMIN || idealEnergyChange >= -LSENERGYTOL) {
            
            if(energyChange < 0) {
                //see if we can return LAMBDAMIN
                if(FFM.computeEnergy(LAMBDAMIN) <= energyCurrent) return LAMBDAMIN;
                else return lambda;
            }
            //can't, just return 0
            else return 0.0;
        }
    }
}

