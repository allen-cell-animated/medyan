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
double BoundaryRepulsion<BRepulsionInteractionType>::ComputeEnergy(BoundaryElement* pbe, double d) {
    double U = 0.0;
    double k_rep = pbe->getRepulsionConst();
    double screenLength = pbe->getScreeningLength();
    
    if (d == 0.0){
        
        for (auto pb: pbe->beads()){
           U+= _FFType.ComputeEnergy(pb, pbe->distance(pb->coordinate), k_rep, screenLength);
        }
    }
    else {
        for (auto pb: pbe->beads()){
            U+=_FFType.ComputeEnergy(pb, pbe->stretchedDistance(pb->coordinate, pb->force, d), k_rep, screenLength);
        }
    }
    return U;
}

template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::ComputeForces(BoundaryElement* pbe) {
    
    double k_rep = pbe->getRepulsionConst();
    double screenLength = pbe->getScreeningLength();
    
    for (auto pb: pbe->beads()){
        auto normal = pbe->normal();
        _FFType.ComputeForces(pb, pbe->distance(pb->coordinate), normal, k_rep, screenLength);
    }
}


template <class BRepulsionInteractionType>
void BoundaryRepulsion<BRepulsionInteractionType>::ComputeForcesAux(BoundaryElement* pbe) { /// Needed for Conjugated Gradient minimization;
    
    double k_rep = pbe->getRepulsionConst();
    double screenLength = pbe->getScreeningLength();
    
    for (auto pb: pbe->beads()){
        auto normal = pbe->normal();
        _FFType.ComputeForcesAux(pb, pbe->distance(pb->coordinate), normal, k_rep, screenLength);
    }
}

///Template specializations
template double BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeEnergy(BoundaryElement* pbe, double d);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeForces(BoundaryElement* pbe);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeForcesAux(BoundaryElement* pbe);
template double BoundaryRepulsion<BoundaryRepulsionExp>::ComputeEnergy(BoundaryElement* pbe, double d);
template void BoundaryRepulsion<BoundaryRepulsionExp>::ComputeForces(BoundaryElement* pbe);
template void BoundaryRepulsion<BoundaryRepulsionExp>::ComputeForcesAux(BoundaryElement* pbe);