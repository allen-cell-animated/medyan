
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

#include "Mechanics/ForceField/Membrane/MembraneFF.h"

#include "Mechanics/ForceField/Membrane/MembraneStretching.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingHarmonic.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingVoronoiHarmonic.h"

#include "Mechanics/ForceField/Membrane/MembraneBending.h"
#include "Mechanics/ForceField/Membrane/MembraneBendingVoronoiHelfrich.h"

#include "Structure/SurfaceMesh/Membrane.hpp"

MembraneFF::MembraneFF (string& stretching, string& bending) {
    
    if (stretching == "TRIANGLE")
        _membraneInteractionVector.emplace_back(
            new MembraneStretching<MembraneStretchingHarmonic>()
        );
    else if(stretching == "VORONOI")
        _membraneInteractionVector.emplace_back(
            new MembraneStretching<MembraneStretchingVoronoiHarmonic>()
        );
    else if(stretching == "") {}
    else {
        cout << "Membrane stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (bending == "HELFRICH")
        _membraneInteractionVector.emplace_back(
            new MembraneBending<MembraneBendingVoronoiHelfrich>()
        );
    else if(bending == "") {}
    else {
        cout << "Membrane bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void MembraneFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit membrane..." << endl;
    _culpritInteraction->_membraneCulprit->printSelf();
    
    cout << endl;
}


double MembraneFF::computeEnergy(double* coord, bool stretched) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _membraneInteractionVector) {
        
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

void MembraneFF::computeForces(double* coord, double* f) {
    
    for (auto &interaction : _membraneInteractionVector)
        interaction->computeForces(coord, f);
}
