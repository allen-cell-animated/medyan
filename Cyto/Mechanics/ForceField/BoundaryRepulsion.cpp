//
//  BoundaryRepulsion.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsion.h"
#include "BoundaryRepulsionLJ.h"
#include "BoundaryRepulsionExp.h"
#include "BoundaryElement.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
double BoundaryRepulsion<BRepulsionInteractionType>::computeEnergy(BoundaryElement* be, double d) {
    double U = 0.0;
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    if (d == 0.0){
        for (auto b: be->getBeads())
           U+= _FFType.computeEnergy(b, abs(be->distance(b->coordinate)), k_rep, screenLength);
    }
    else {
        for (auto b: be->getBeads())
            U+=_FFType.computeEnergy(b, abs(be->stretchedDistance(b->coordinate, b->force, d)), k_rep, screenLength);
    }
    
    return U;
}

template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::computeForces(BoundaryElement* be) {
    
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    for (auto b: be->getBeads()){
        auto normal = be->normal(b->coordinate);
        _FFType.computeForces(b, abs(be->distance(b->coordinate)), normal, k_rep, screenLength);
    }
}


template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::computeForcesAux(BoundaryElement* be) { /// Needed for Conjugated Gradient minimization;
    
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    for (auto b: be->getBeads()){
        auto normal = be->normal(b->coordinateAux);
        _FFType.computeForcesAux(b, abs(be->distance(b->coordinateAux)), normal, k_rep, screenLength);
    }
}

///Template specializations
template double BoundaryRepulsion<BoundaryRepulsionLJ>::computeEnergy(BoundaryElement* be, double d);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::computeForces(BoundaryElement* be);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::computeForcesAux(BoundaryElement* be);
template double BoundaryRepulsion<BoundaryRepulsionExp>::computeEnergy(BoundaryElement* be, double d);
template void BoundaryRepulsion<BoundaryRepulsionExp>::computeForces(BoundaryElement* be);
template void BoundaryRepulsion<BoundaryRepulsionExp>::computeForcesAux(BoundaryElement* be);