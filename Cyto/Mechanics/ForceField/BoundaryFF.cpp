//
//  BoundaryFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryFF.h"

#include "BoundaryRepulsion.h"
#include "BoundaryRepulsionLJ.h"
#include "BoundaryRepulsionExp.h"

#include "BoundaryElementDB.h"
#include "Bead.h"

BoundaryFF::BoundaryFF (string interaction1, string interaction2, string interaction3) {
    
    if (interaction1 == "REPULSIONLJ") {_BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionLJ>());}
    if (interaction1 == "REPULSIONEXP") {_BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionExp>());}
}


double BoundaryFF::computeEnergy(double d) {
    double U = 0;
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance()) {
            
            for(auto &neighbor : neighborList->getNeighbors(boundaryElement)) {
                Bead* bead = static_cast<Bead*>(neighbor);
                U += boundaryInteraction->computeEnergy(boundaryElement, bead, d);
            }
        }
    }
    return U;
}

void BoundaryFF::computeForces() {
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance()) {
            
            for(auto &neighbor : neighborList->getNeighbors(boundaryElement)) {
                Bead* bead = static_cast<Bead*>(neighbor);
                boundaryInteraction->computeForces(boundaryElement, bead);
            }
        }
    }
}

void BoundaryFF::computeForcesAux() {
    
    for (auto &boundaryInteraction : _BoundaryInteractionVector){
        
        auto neighborList = boundaryInteraction->getNeighborList();
        for (auto boundaryElement: *BoundaryElementDB::instance()) {
            
            for(auto &neighbor : neighborList->getNeighbors(boundaryElement)) {
                Bead* bead = static_cast<Bead*>(neighbor);
                boundaryInteraction->computeForcesAux(boundaryElement, bead);
            }
        }
    }
}