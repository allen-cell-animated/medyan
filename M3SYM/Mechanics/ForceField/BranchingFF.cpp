//
//  BranchingFF.cpp
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BranchingFF.h"

#include "BranchingPointDB.h"

#include "BranchingStretching.h"
#include "BranchingStretchingHarmonic.h"

#include "BranchingBending.h"
#include "BranchingBendingCosine.h"

#include "BranchingDihedral.h"
#include "BranchingDihedralCosine.h"

BranchingFF::BranchingFF (string& stretching, string& bending, string& dihedral)
{
    if (stretching == "HARMONIC")
        _branchingInteractionVector.emplace_back(new BranchingStretching<BranchingStretchingHarmonic>());
    if (bending == "COSINE")
        _branchingInteractionVector.emplace_back(new BranchingBending<BranchingBendingCosine>());
//  if (dihedral == "COSINE")
//      _branchingInteractionVector.emplace_back(new BranchingDihedral<BranchingDihedralCosine>());
//
//    if (Twisting == "HARMONIC")
//        _filamentInteractionVector.push_back(new FilamentTwisting<FilamentTwistingHarmonic>());
}


double BranchingFF::computeEnergy(double d) {
    double U_branch = 0;
    for (auto branch: *BranchingPointDB::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            U_branch += branchingInteraction.get()->computeEnergy(branch, d);
    return U_branch;
}

void BranchingFF::computeForces() {
    for (auto branch: *BranchingPointDB::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            branchingInteraction.get()->computeForces(branch);
}

void BranchingFF::computeForcesAux() {
    
    for (auto branch: *BranchingPointDB::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            branchingInteraction.get()->computeForcesAux(branch);
}
