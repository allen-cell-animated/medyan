
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

#include "ForceFieldManager.h"

double ForceFieldManager::computeEnergy(double d) {
    
    double energy = 0;
    for(auto &f : _forceFields)
        energy += f->computeEnergy(d);
    
    //if energy is infinity, exit ungracefully.
    if(energy == numeric_limits<double>::infinity()) {
        cout <<
        "Energy became infinite. Try adjusting equilibration step size." << endl;
        exit(EXIT_FAILURE);
    }
    return energy;
}

void ForceFieldManager::computeForces() {
    resetForces();
    for(auto &f : _forceFields) f->computeForces();
}

void ForceFieldManager::computeForcesAux() {
    resetForcesAux();
    for(auto &f : _forceFields)
        f->computeForcesAux();
}

void ForceFieldManager::resetForces() {
    
    for(auto it: *BeadDB::instance())
        it->force.assign (3, 0); //Set force to zero;
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto it: *BeadDB::instance())
        it->forceAux.assign (3, 0); //Set forceAux to zero;
}