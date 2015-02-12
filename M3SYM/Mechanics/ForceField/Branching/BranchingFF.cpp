
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
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
    else if(stretching == "") {}
    else {
        cout << "Branching stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(bending == "COSINE")
        _branchingInteractionVector.emplace_back(
        new BranchingBending<BranchingBendingCosine>());
    else if(bending == "") {}
    else {
        cout << "Branching bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(dihedral == "COSINE")
      _branchingInteractionVector.emplace_back(
      new BranchingDihedral<BranchingDihedralCosine>());
    else if(dihedral == "") {}
    else {
        cout << "Branching dihedral FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(position == "COSINE")
      _branchingInteractionVector.emplace_back(new
          BranchingPosition<BranchingPositionCosine>());
    else if(position == "") {}
    else {
        cout << "Branching position FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }

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
