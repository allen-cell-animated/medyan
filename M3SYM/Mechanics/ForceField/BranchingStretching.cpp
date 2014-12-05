//
//  BranchingStretching.cpp
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BranchingStretching.h"


#include "BranchingStretchingHarmonic.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BStretchingInteractionType>
double BranchingStretching<BStretchingInteractionType>::computeEnergy(BranchingPoint* b, double d) {
    
//    Bead* b1 =
//    Bead* b2 =
//    Bead* b3 =
//    double kStretch =
//    double l0 =
//    
//    double position =
//
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, position, kStretch, l0);
    else
        return _FFType.energy(b1, b2, b3, position, kStretch, l0, d);
    
}

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForces(BranchingPoint* b) {
    
    //    Bead* b1 =
    //    Bead* b2 =
    //    Bead* b3 =
    //    double kStretch =
    //    double l0 =
    //
    //    double position =
    //
    
    _FFType.forces(b1, b2, b3, position, kStretch, l0);
    
}


template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForcesAux(BranchingPoint* b) {
    
    //    Bead* b1 =
    //    Bead* b2 =
    //    Bead* b3 =
    //    double kStretch =
    //    double l0 =
    //
    //    double position =
    //
    
    _FFType.forcesAux(b1, b2, b3, position, kStretch, l0);
    
}


///Template specializations
template double BranchingStretching<BranchingStretchingHarmonic>::computeEnerg(BranchingPoint* b, double d);
template void  BranchingStretching<BranchingStretchingHarmonic>::computeForces(BranchingPoint* b);
template void  BranchingStretching<BranchingStretchingHarmonic>::computeForcesAux(BranchingPoint* b);