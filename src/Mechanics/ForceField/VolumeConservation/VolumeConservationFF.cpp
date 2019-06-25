
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

#include "Mechanics/ForceField/VolumeConservation/VolumeConservationFF.h"

#include "Mechanics/ForceField/VolumeConservation/VolumeConservationMembrane.h"
#include "Mechanics/ForceField/VolumeConservation/VolumeConservationMembraneHarmonic.h"

#include "Structure/SurfaceMesh/Membrane.hpp"

VolumeConservationFF::VolumeConservationFF(string& type) {
    if (type == "MEMBRANE")
        _volumeConservationInteractionVector.emplace_back(
            new VolumeConservationMembrane<VolumeConservationMembraneHarmonic>());
    else if(type == "") {}
    else {
        cout << "Volume conservation FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void VolumeConservationFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit membrane..." << endl;
    _culpritInteraction->_membraneCulprit->printSelf();
    
    cout << endl;
}

floatingpoint VolumeConservationFF::computeEnergy(floatingpoint* coord, bool stretched) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _volumeConservationInteractionVector) {
        
        U_i = interaction->computeEnergy(coord, stretched);
                
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;
        
    }
    return U;
}

void VolumeConservationFF::computeForces(floatingpoint* coord, floatingpoint* f) {
    
    for (auto &interaction : _volumeConservationInteractionVector)
        interaction->computeForces(coord, f);
}
