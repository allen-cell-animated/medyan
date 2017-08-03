
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

#include "ForceFieldManager.h"

#include "CGMethod.h"

double ForceFieldManager::computeEnergy(double *coord, double *f, double d, bool verbose) {
    
    double energy = 0;
    for(auto &ff : _forceFields) {
        
        auto tempEnergy = ff->computeEnergy(coord, f, d);
        
        if(verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;
        
        //if energy is infinity, exit with infinity.
        if(tempEnergy <= -1) {
            
            //if this is the current energy, exit ungracefully
            if(d == 0.0) {
                
                cout << "Energy = " << tempEnergy << endl;
                
                cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
                cout << "The culprit was ... " << ff->getName() << endl;
                
                //get the culprit in output
                ff->whoIsCulprit();
                
                exit(EXIT_FAILURE);
            }
            //if this is a minimization try, just return infinity
            else return numeric_limits<double>::infinity();
        }
        else energy += tempEnergy;
        
    }
    return energy;
}

void ForceFieldManager::computeForces(double *coord, double *f) {
    
    //reset
    resetForces();
    
    //recompute
    for(auto &ff : _forceFields) ff->computeForces(coord, f);

    //WILL HAVE TO COPY AUXS AFTER THIS CALL
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
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto b: Bead::getBeads())
        b->forceAux.assign (3, 0); //Set forceAux to zero;
}

void ForceFieldManager::copyForces(double *fprev, double *f) {
    
    for (int i = 0; i < CGMethod::N; i++)
        fprev[i] = f[i];
}
