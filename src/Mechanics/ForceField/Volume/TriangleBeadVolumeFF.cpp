
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

#include "Mechanics/ForceField/Volume/TriangleBeadVolumeFF.hpp"

#include "Mechanics/ForceField/Volume/TriangleBeadExclVolume.hpp"
#include "Mechanics/ForceField/Volume/TriangleCylinderBeadExclVolRepulsion.hpp"
#include "Structure/Cylinder.h"
#include "Structure/SurfaceMesh/Triangle.hpp"

TriangleBeadVolumeFF::TriangleBeadVolumeFF (string& type) {
    if (type == "REPULSION")
        _triangleBeadVolInteractionVector.emplace_back(
        new TriangleBeadExclVolume <TriangleCylinderBeadExclVolRepulsion>());
    else if(type == "") {}
    else {
        cout << "Volume FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void TriangleBeadVolumeFF::vectorize() {
    for(auto& interaction : _triangleBeadVolInteractionVector)
        interaction->vectorize();
}

void TriangleBeadVolumeFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit triangle and bead..." << endl;
    _culpritInteraction->triangleCulprit_->printSelf();
    _culpritInteraction->beadCulprit_    ->printSelf();
    
    cout << endl;
}

floatingpoint TriangleBeadVolumeFF::computeEnergy(floatingpoint* coord, bool stretched) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _triangleBeadVolInteractionVector) {
        
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

void TriangleBeadVolumeFF::computeForces(floatingpoint* coord, floatingpoint* f) {
    
    for (auto &interaction : _triangleBeadVolInteractionVector)
        interaction->computeForces(coord, f);
}

void TriangleBeadVolumeFF::computeLoadForces() {
    for(auto& interaction: _triangleBeadVolInteractionVector)
        interaction->computeLoadForces();
}
void TriangleBeadVolumeFF::computeLoadForce(Cylinder* c, LoadForceEnd end) const {
    for(const auto& interaction : _triangleBeadVolInteractionVector)
        interaction->computeLoadForce(c, end);
}

vector<NeighborList*> TriangleBeadVolumeFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _triangleBeadVolInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}

