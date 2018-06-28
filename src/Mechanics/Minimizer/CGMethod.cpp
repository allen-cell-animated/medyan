
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

Bead* CGMethod::maxBead() {
    
    double maxF = 0;
    double force;
    Bead* retbead = nullptr;
    
    //calc max force
    for(auto b: Bead::getBeads()) {
        
        force = b->FADotFA();
        
        if(force >= maxF) {
            retbead = b;
            maxF = force;
        }
    }
    return retbead;
}

void CGMethod::startMinimization() {

    for(auto b: Bead::getBeads())
        b->coordinateB = b->coordinateP = b->coordinate;
}


void CGMethod::moveBeads(double d)
{
	for(auto b: Bead::getBeads()) {
        if(!b->getstaticstate()){
        b->coordinate[0] = b->coordinate[0] + d * b->force[0];
        b->coordinate[1] = b->coordinate[1] + d * b->force[1];
        b->coordinate[2] = b->coordinate[2] + d * b->force[2];
        }

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
    int iter = 0;
    //return zero if no forces
    if(f == 0.0){std::cout<<"lambda determined in "<<iter<<endl;return 0.0;}
    
    //calculate first lambda
    double lambda = min(LAMBDAMAX, MAXDIST / f);
    double currentEnergy = FFM.computeEnergy(0.0);
    
    //backtracking loop
    while(true) {
        iter++;
        
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        
        double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA();
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= idealEnergyChange) {
            std::cout<<"SL2 BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA "<<allFDotFA()<<endl;
            std::cout<<"SL2 energyChange "<<energyChange<<" idealEnergyChange "<<idealEnergyChange<<" lambda "<<lambda<<endl;
            std::cout<<"lambda determined in "<<iter<<endl; return lambda;}
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        if(lambda <= 0.0 || lambda <= LAMBDATOL)
        {
            std::cout<<"SL2 BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA "<<allFDotFA()<<endl;
            std::cout<<"SL2 energyChange "<<energyChange<<"idealEnergyChange "<<idealEnergyChange<<" lambda "<<lambda<<endl;
            std::cout<<"lambda determined in "<<iter<<endl;
            return 0.0;}
        std::cout<<"SL2 BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA "<<allFDotFA()<<endl;
        std::cout<<"SL2 energyChange "<<energyChange<<"idealEnergyChange "<<idealEnergyChange<<" lambda "<<lambda<<endl;
    }
}

double CGMethod::safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                                    double LAMBDAMAX) {
    
    //reset safe mode
    _safeMode = false;
    
    //calculate first lambda
    double lambda = LAMBDAMAX;
    double currentEnergy = FFM.computeEnergy(0.0);
    int iter = 0;
    //backtracking loop
    while(true) {
        iter++;
        std::cout<<"safe z"<<endl;
        //new energy when moved by lambda
        double energyLambda = FFM.computeEnergy(lambda);
        double energyChange = energyLambda - currentEnergy;
        
        //return if ok
        if(energyChange <= 0.0) {
            std::cout<<"safe energyChange "<<energyChange<<" lambda "<<lambda<<endl;
            std::cout<<"lambda determined in "<<iter<<endl;
            return lambda;}
        
        //reduce lambda
        lambda *= LAMBDAREDUCE;
        
        //just shake if we cant find an energy min,
        //so we dont get stuck
        if(lambda <= 0.0 || lambda <= LAMBDATOL){
            std::cout<<"safe energyChange "<<energyChange<<" lambda "<<lambda<<endl;
            std::cout<<"lambda determined in "<<iter<<endl;
            return MAXDIST / maxF();
        }
    }
}
