
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
}

void ForceFieldManager::computeForcesAux() {
    
    //reset just aux
    resetForcesAux();
    
    //recompute
    for(auto &f : _forceFields) f->computeForcesAux();
}

void ForceFieldManager::computeForcesAuxP() {
    
    //copy to auxp
    for(auto b: Bead::getBeads())
        b->forceAuxP = b->forceAux;
}

void ForceFieldManager::resetForces() {
    
    for(auto b: Bead::getBeads())
        b->force.assign (3, 0); //Set force to zero;
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto b: Bead::getBeads())
        b->forceAux.assign (3, 0); //Set forceAux to zero;
}