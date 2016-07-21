
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ForceFieldManager.h"

#include "BoundaryElement.h"

double ForceFieldManager::computeEnergy(double d, bool verbose) {
    
    double energy = 0;
    for(auto &f : _forceFields) {
        
        auto tempEnergy = f->computeEnergy(d);
        
        if(verbose) cout << f->getName() << " energy = " << tempEnergy << endl;
        
        //if energy is infinity, exit with infinity.
        if(tempEnergy <= -1) {
            
            //if this is the current energy, exit ungracefully
            if(d == 0.0) {
                
                cout << "Energy = " << tempEnergy << endl;
                
                cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
                cout << "The culprit was ... " << f->getName() << endl;
                
                //get the culprit in output
                f->whoIsCulprit();
                
                exit(EXIT_FAILURE);
            }
            //if this is a minimization try, just return infinity
            else return numeric_limits<double>::infinity();
        }
        else energy += tempEnergy;
        
    }
    return energy;
}

void ForceFieldManager::computeForces() {
    
    //reset
    resetForces();
    
    //recompute
    for(auto &f : _forceFields) f->computeForces();
    
    //copy to auxs
    for(auto b: Bead::getBeads())
        b->forceAux = b->forceAuxP = b->force;

    for(auto be: BoundaryElement::getBoundaryElements())
    	be->boundary_forceAux = be->boundary_forceAuxP = be->boundary_force;

}

void ForceFieldManager::computeForcesAux() {
    
    //reset just aux
    resetForcesAux();
    
    //recompute
    for(auto &f : _forceFields)
        f->computeForcesAux();

}

void ForceFieldManager::computeForcesAuxP() {
    
    //copy to auxp
    for(auto b: Bead::getBeads())
        b->forceAuxP = b->forceAux;

    //added by jl135, for storing previous boundary_forceAux
    for(auto be: BoundaryElement::getBoundaryElements())
    	be->boundary_forceAuxP = be->boundary_forceAux;

}

void ForceFieldManager::computeLoadForces() {
    
    for(auto &f : _forceFields)
        f->computeLoadForces();
    
    //reset lfi as well
    for(auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}


void ForceFieldManager::resetForces() {
    
    for(auto b: Bead::getBeads()) {
        b->force.assign (3, 0); //Set force to zero;
        std::memset((void*)(&b->loadForcesP[0]), 0, sizeof(b->loadForcesP));  //Set load force to zero;
        std::memset((void*)(&b->loadForcesM[0]), 0, sizeof(b->loadForcesM));  //Set load force to zero;
    }

    //set boundary force to zero, added by jl135
    for(auto be: BoundaryElement::getBoundaryElements())
    	be->boundary_force = 0;
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto b: Bead::getBeads())
        b->forceAux.assign (3, 0); //Set forceAux to zero;

    //added by jl135, reset force field for boundary
    for(auto be: BoundaryElement::getBoundaryElements())
    	be->boundary_forceAux = 0;
}



