
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

#include "BranchingFF.h"

#include "BranchingStretching.h"
#include "BranchingStretchingHarmonic.h"

#include "BranchingBending.h"
#include "BranchingBendingCosine.h"

#include "BranchingDihedral.h"
#include "BranchingDihedralCosine.h"

#include "BranchingPosition.h"
#include "BranchingPositionCosine.h"

#include "BranchingPoint.h"

BranchingFF::BranchingFF(string& stretching, string& bending,
                         string& dihedral, string& position)
{
    if(stretching == "HARMONIC")
        _branchingInteractionVector.emplace_back(
        new BranchingStretching<BranchingStretchingHarmonic>());
    
    if(bending == "COSINE")
        _branchingInteractionVector.emplace_back(
        new BranchingBending<BranchingBendingCosine>());
    
    if(dihedral == "COSINE")
      _branchingInteractionVector.emplace_back(
      new BranchingDihedral<BranchingDihedralCosine>());
    
    if(position == "COSINE")
      _branchingInteractionVector.emplace_back(new
          BranchingPosition<BranchingPositionCosine>());

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
