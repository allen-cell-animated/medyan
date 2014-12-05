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
#include "BranchingBendingHarmonic.h"

#include "BranchingDihedral.h"
#include "BranchingDihedralHarmonic.h"

BranchingFF::BranchingFF (string& stretching, string& bending, string& dihedral)
{
    if (stretching == "HARMONIC")
        _branchingInteractionVector.emplace_back(new BranchingStretching<BranchingStretchingHarmonic>());
    if (bending == "HARMONIC")
        _branchingInteractionVector.emplace_back(new BranchingBending<BranchingBendingHarmonic>());
    if (bending == "HARMONIC")
        _branchingInteractionVector.emplace_back(new BranchingDihedral<BranchingDihedralHarmonic>());
    
    
    //if (Twisting == "HARMONIC") {_filamentInteractionVector.push_back(new FilamentTwisting<FilamentTwistingHarmonic>());}
}


double BranchingFF::computeEnergy(double d) {
    double U_branch = 0;
    for (auto branch: *BranchingPointDB::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            U_branch += branchingInteraction.get()->computeEnergy(branch, d);
    return U_branch;
}

void BranchingFF::computeForces() {
    for (auto branch: *BranchingPoint::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            branchingInteraction.get()->computeForces(fil);
}

void FilamentFF::computeForcesAux() {
    
    for (auto branch: *BranchingPoint::instance())
        for (auto &branchingInteraction : _branchingInteractionVector)
            branchingInteraction.get()->computeForcesAux(fil);
}
