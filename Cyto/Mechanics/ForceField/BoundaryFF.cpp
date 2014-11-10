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

BoundaryFF::BoundaryFF (string interaction1, string interaction2, string interaction3) {
    
    if (interaction1 == "REPULSIONLJ") {_BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionLJ>());}
    if (interaction1 == "REPULSIONEXP") {_BoundaryInteractionVector.emplace_back(new BoundaryRepulsion<BoundaryRepulsionExp>());}
}


double BoundaryFF::computeEnergy(double d) {
    double U_bound = 0;
    
    for ( auto it: *BoundaryElementDB::instance(BoundaryElementDBKey()) ) {
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
            U_bound += boundaryInteraction.get()->computeEnergy(it, d);
    }
    return U_bound;
}

void BoundaryFF::computeForces() {
    for ( auto it: *BoundaryElementDB::instance(BoundaryElementDBKey()) ) {
        
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
            boundaryInteraction.get()->computeForces(it);
    }
}

void BoundaryFF::computeForcesAux() {
    for ( auto it: *BoundaryElementDB::instance(BoundaryElementDBKey()) ) {
        
        for (auto &boundaryInteraction : _BoundaryInteractionVector)
           boundaryInteraction.get()->computeForcesAux(it);
    }
}