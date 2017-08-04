
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
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit motor..." << endl;
    _culpritInteraction->_motorCulprit->printSelf();
    
    cout << endl;
}

void MotorGhostFF::vectorize() {
    
    for (auto &interaction : _motorGhostInteractionVector)
        interaction->vectorize();
}

void MotorGhostFF::cleanup() {
    
    for (auto &interaction : _motorGhostInteractionVector)
        interaction->deallocate();
}


double MotorGhostFF::computeEnergy(double *coord, double *f, double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _motorGhostInteractionVector) {
        
        U_i = interaction->computeEnergy(coord, f, d);
        
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;
        
    }
    return U;
}

void MotorGhostFF::computeForces(double *coord, double *f) {
    
    for (auto &interaction : _motorGhostInteractionVector)
        interaction->computeForces(coord, f);
}

