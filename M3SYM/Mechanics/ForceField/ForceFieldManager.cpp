
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

#include "ForceFieldManager.h"

double ForceFieldManager::computeEnergy(double d) {
    
    double energy = 0;
    for(auto &f : _forceFields) {
        auto tempEnergy = f->computeEnergy(d);
        
        //if energy is infinity, exit with infinity.
        if(tempEnergy == -1) {
            cout << "Warning: Energy became garbage. Try adjusting equilibration step size."
                 << endl;
            cout << "The culprit was... " << f->getName() << endl;
            return numeric_limits<double>::infinity();
        }
        else {
            energy += tempEnergy;
        }
    }
    return energy;
}

void ForceFieldManager::computeForces() {
    
    //reset
    resetForces();
    
    //recompute
    for(auto &f : _forceFields) f->computeForces();
    
    //copy to aux
    for(auto it: *BeadDB::instance())
        it->forceAux = it->force;
}

void ForceFieldManager::computeForcesAux() {
    
    //reset just aux
    resetForcesAux();
    
    //recompute
    for(auto &f : _forceFields) f->computeForcesAux();
}

void ForceFieldManager::resetForces() {
    
    for(auto it: *BeadDB::instance())
        it->force.assign (3, 0); //Set force to zero;
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto it: *BeadDB::instance())
        it->forceAux.assign (3, 0); //Set forceAux to zero;
}