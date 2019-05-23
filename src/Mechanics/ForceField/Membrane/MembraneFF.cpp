
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

#include <stdexcept>

#include "Mechanics/ForceField/Membrane/MembraneStretching.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"

#include "Mechanics/ForceField/Membrane/MembraneBending.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "util/io/log.h"

MembraneFF::MembraneFF (const string& stretching, const string& stretchingAccu, const string& bending) {
    
    if (stretching == "HARMONIC")
        if(stretchingAccu == "TRIANGLE")
            _membraneInteractionVector.emplace_back(
                new MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByTriangle >()
            );
        else if(stretchingAccu == "VERTEX")
            _membraneInteractionVector.emplace_back(
                new MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByVertex >()
            );
        else {
            LOG(ERROR) << "Membrane stretching accumulation type " << stretchingAccu << " is not recognized.";
            throw std::runtime_error("Membrane stretching accumulation type not recognized");
        }

    else if(stretching == "LINEAR")
        if(stretchingAccu == "TRIANGLE")
            _membraneInteractionVector.emplace_back(
                new MembraneStretching< MembraneStretchingLinear, MembraneStretchingAccumulationType::ByTriangle >()
            );
        else if(stretchingAccu == "VERTEX")
            _membraneInteractionVector.emplace_back(
                new MembraneStretching< MembraneStretchingLinear, MembraneStretchingAccumulationType::ByVertex >()
            );
        else {
            LOG(ERROR) << "Membrane stretching accumulation type " << stretchingAccu << " is not recognized.";
            throw std::runtime_error("Membrane stretching accumulation type not recognized");
        }

    else if(stretching == "") {}

    else {
        LOG(ERROR) << "Membrane stretching FF type " << stretching << " is not recognized.";
        throw std::runtime_error("Membrane stretching FF type not recognized");
    }
    
    if (bending == "HELFRICH")
        _membraneInteractionVector.emplace_back(
            new MembraneBending<MembraneBendingHelfrich>()
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
