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
#include "BoundaryElementDB.h"

BoundaryFF::BoundaryFF (std::string Interaction1, std::string Interaction2, std::string Interaction3) {
    
    if (Interaction1 == "REPULSIONLJ") {_BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionLJ>());}
}


double BoundaryFF::ComputeEnergy(double d) {
    double U_bound = 0;
    
    for ( auto it: *BoundaryElementDB::Instance(BoundaryElementDBKey()) ) {
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
            U_bound += boundaryInteraction.get()->ComputeEnergy(it, d);
    }
    return U_bound;
}

void BoundaryFF::ComputeForces() {
    for ( auto it: *BoundaryElementDB::Instance(BoundaryElementDBKey()) ) {
        
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
            boundaryInteraction.get()->ComputeForces(it);
    }
}

void BoundaryFF::ComputeForcesAux() {
    for ( auto it: *BoundaryElementDB::Instance(BoundaryElementDBKey()) ) {
        
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
           boundaryInteraction.get()->ComputeForcesAux(it);
    }
}