
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
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    double energy = 0;
    for(auto &f : _forceFields) {
        
        auto tempEnergy = f->computeEnergy(d);
        
        //if energy is infinity, exit with infinity.
        if(tempEnergy == -1) {
            
            //if this is the current energy, exit ungracefully
            if(d == 0.0) {
                cout << "Energy became infinite. Try adjusting minimization parameters." << endl;
                cout << "The culprit was ... " << f->getName() << endl;
                exit(EXIT_FAILURE);
            }
            //if this is a minimization try, just return infinity
            else
                return numeric_limits<double>::infinity();
        }
        else energy += tempEnergy;
    }
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "Energy calculation took " << elapsed_run.count() << endl;
    return energy;
}

void ForceFieldManager::computeForces() {
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    //reset
    resetForces();
    
    //recompute
    for(auto &f : _forceFields) f->computeForces();
    
    //copy to auxs
    for(auto b: *BeadDB::instance())
        b->forceAux = b->forceAuxP = b->force;
    
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "Force calculation took " << elapsed_run.count() << endl;
}

void ForceFieldManager::computeForcesAux() {
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    //reset just aux
    resetForcesAux();
    
    //recompute
    for(auto &f : _forceFields) f->computeForcesAux();
    
    
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "ForceAux calculation took " << elapsed_run.count() << endl;
}

void ForceFieldManager::computeForcesAuxP() {
    
    //copy to auxp
    for(auto b: *BeadDB::instance())
        b->forceAuxP = b->forceAux;
}

void ForceFieldManager::resetForces() {
    
    for(auto it: *BeadDB::instance())
        it->force.assign (3, 0); //Set force to zero;
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto it: *BeadDB::instance())
        it->forceAux.assign (3, 0); //Set forceAux to zero;
}