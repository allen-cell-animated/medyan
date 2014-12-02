
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundaryFF.h"

#include "BoundaryRepulsion.h"
#include "BoundaryRepulsionLJ.h"
#include "BoundaryRepulsionExp.h"

#include "BoundaryElement.h"
#include "Bead.h"

BoundaryFF::BoundaryFF (string interaction1, string interaction2, string interaction3) {
    
    if (interaction1 == "REPULSIONLJ")
        _BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionLJ>());
    if (interaction1 == "REPULSIONEXP")
        _BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionExp>());
}

double BoundaryFF::computeEnergy(double d) {
    double U = 0;
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance())
            for(auto &bead : neighborList->getNeighbors(boundaryElement))
                U += boundaryInteraction->computeEnergy(boundaryElement, bead, d);
    }
    return U;
}

void BoundaryFF::computeForces() {
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance())
            for(auto &bead : neighborList->getNeighbors(boundaryElement))
                boundaryInteraction->computeForces(boundaryElement, bead);
    }
}

void BoundaryFF::computeForcesAux() {
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance())
            for(auto &bead : neighborList->getNeighbors(boundaryElement))
                boundaryInteraction->computeForcesAux(boundaryElement, bead);
    }
}
