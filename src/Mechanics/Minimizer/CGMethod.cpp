
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CGMethod.h"

#include "ForceFieldManager.h"
#include "Bead.h"

long CGMethod::N = 0;

double CGMethod::allFDotF()
{
    double g = 0;
    for(int i = 0; i < N; i++)
        g += force[i] * force[i];
    return g;
}

double CGMethod::allFADotFA()
{
    double g = 0;
    for(int i = 0; i < N; i++)
        g += forceAux[i] * forceAux[i];
    return g;
}

double CGMethod::allFADotFAP()
{
    double g = 0;
    for(int i = 0; i < N; i++)
        g += forceAux[i] * forceAuxPrev[i];
    return g;
}

double CGMethod::allFDotFA()
{
    double g = 0;
    for(int i = 0; i < N; i++)
        g += force[i] * forceAux[i];
    return g;
}

double CGMethod::maxF() {
    
    double maxF = 0;
    double mag;
    
    for(int i = 0; i < N; i++) {
        mag = sqrt(forceAux[i]*forceAux[i]);
        if(mag > maxF) maxF = mag;
    }
    
    return maxF;
}

Bead* CGMethod::maxBead() {
    
    double maxF = 0;
    double currentF;
    long index = 0;
    
    for (int i = 0; i < N; i++) {
        
        currentF = forceAux[i] * forceAux[i];
        if(currentF > maxF) {
            index = i;
            maxF = currentF;
        }
    }
    return Bead::getBeads()[index];
}

void CGMethod::moveBeads(double d)
{
    ///<NOTE: Ignores static beads for now.
    //if(!b->getstaticstate())
    
    for (int i = 0; i < N; i++) {
        coord[i] = coord[i] + d * force[i];
    }
    
}

void CGMethod::shiftGradient(double d)
{
    for (int i = 0; i < N; i ++)
        force[i] = forceAux[i] + d * force[i];
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

void CGMethod::startMinimization() {
 
    //COPY BEAD DATA
    N = 3 * Bead::getBeads().size();
    allocate(N);

    //coord management
    long i = 0;
    long index = 0;
    for(auto b: Bead::getBeads()) {
        
        //set bead index
        b->_dbIndex = i;
        
        //flatten indices
        index = 3 * i;
        coord[index] = b->coordinate[0];
        coord[index + 1] = b->coordinate[1];
        coord[index + 2] = b->coordinate[2];
        
        b->coordinateP = b->coordinate;
        i++;
    }
}

void CGMethod::endMinimization() {
    
    ///RECOPY BEAD DATA
    //coord management
    long i = 0;
    long index = 0;
    for(auto b: Bead::getBeads()) {
        
        //flatten indices
        index = 3 * i;
        b->coordinate[0] = coord[index];
        b->coordinate[1] = coord[index + 1];
        b->coordinate[2] = coord[index + 2];
        
        i++;
    }
    
    deallocate();
}


double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                double LAMBDAMAX) {

    double f = maxF();
    
    //return zero if no forces
    if(f == 0.0) return 0.0;
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / f);
    double currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(coord, force, lambda);
        
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA();
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || lambda <= LAMBDATOL)
            return 0.0;
    }
}

double CGMethod::safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                    double LAMBDAMAX) {
    
    //reset safe mode
    _safeMode = false;
    
    //calculate first lambda
    double lambda = LAMBDAMAX;
    double currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    
    //backtracking loop
    while(true) {
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(coord, force, lambda);
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= 0.0) return lambda;
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        //just shake if we cant find an energy min,
        //so we dont get stuck
        if(lambda <= 0.0 || lambda <= LAMBDATOL)
            return MAXDIST / maxF();
    }
}
