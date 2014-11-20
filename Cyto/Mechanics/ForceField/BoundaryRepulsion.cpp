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
double BoundaryRepulsion<BRepulsionInteractionType>::computeEnergy(BoundaryElement* be, Bead* b, double d) {
    double U = 0.0;
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    if (d == 0.0)
        U+= _FFType.computeEnergy(b, be->distance(b->coordinate), k_rep, screenLength);
    else
        U+=_FFType.computeEnergy(b, be->stretchedDistance(b->coordinate, b->force, d), k_rep, screenLength);
    
    return U;
}

template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::computeForces(BoundaryElement* be, Bead* b) {
    
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    auto normal = be->normal(b->coordinate);
    _FFType.computeForces(b, be->distance(b->coordinate), normal, k_rep, screenLength);
    
}


template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::computeForcesAux(BoundaryElement* be, Bead* b) {
    
    double k_rep = be->getRepulsionConst();
    double screenLength = be->getScreeningLength();
    
    auto normal = be->normal(b->coordinateAux);
    _FFType.computeForcesAux(b, be->distance(b->coordinateAux), normal, k_rep, screenLength);
}

///Template specializations
template double BoundaryRepulsion<BoundaryRepulsionLJ>::computeEnergy(BoundaryElement* be, Bead* b, double d);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::computeForces(BoundaryElement* be, Bead* b);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::computeForcesAux(BoundaryElement* be, Bead* b);
template double BoundaryRepulsion<BoundaryRepulsionExp>::computeEnergy(BoundaryElement* be, Bead* b, double d);
template void BoundaryRepulsion<BoundaryRepulsionExp>::computeForces(BoundaryElement* be, Bead* b);
template void BoundaryRepulsion<BoundaryRepulsionExp>::computeForcesAux(BoundaryElement* be, Bead* b);