
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

#include "MotorGhostFF.h"

#include "MotorGhostStretching.h"
#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"

MotorGhostFF::MotorGhostFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC")
        _motorGhostInteractionVector.emplace_back(
            new MotorGhostStretching<MotorGhostStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Motor stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void MotorGhostFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Printing the culprit motor..." << endl;
    
    _motorCulprit->printInfo();
    
    cout << endl;
}

double MotorGhostFF::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto &interaction : _motorGhostInteractionVector) {
        
        for (auto m: MotorGhost::getMotorGhosts()) {
            
            U_i = interaction->computeEnergy(m, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i) {
                
                //set culprit and return
                _motorCulprit = m;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    return U;
}

void MotorGhostFF::computeForces() {
    
    for (auto &interaction : _motorGhostInteractionVector) {
        
        for (auto m: MotorGhost::getMotorGhosts())
            interaction->computeForces(m);
    }
}

void MotorGhostFF::computeForcesAux() {
    
    for (auto &interaction : _motorGhostInteractionVector) {
        
        for (auto m: MotorGhost::getMotorGhosts())
            interaction->computeForcesAux(m);
    }
}

