
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

#include "MController.h"

#include "SubSystem.h"

#include "FilamentFF.h"
#include "LinkerFF.h"
#include "MotorGhostFF.h"
#include "BoundaryFF.h"
#include "VolumeCylindricalFF.h"
#include "BranchingFF.h"

#include "ConjugateGradient.h"

void MController::initializeMinAlgorithms (MechanicsAlgorithm& MAlgorithm) {
    
    if (MAlgorithm.ConjugateGradient == "FLETCHERRIEVES")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<FletcherRieves>(MAlgorithm.gradientTolerance,
                                              MAlgorithm.maxDistance));
    else if (MAlgorithm.ConjugateGradient == "POLAKRIBIERE")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<PolakRibiere>(MAlgorithm.gradientTolerance,
                                            MAlgorithm.maxDistance));
    else if (MAlgorithm.ConjugateGradient == "STEEPESTDESCENT")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<SteepestDescent>(MAlgorithm.gradientTolerance,
                                               MAlgorithm.maxDistance));
    
    else {
        cout << "Conjugate gradient method not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void MController::initializeFF (MechanicsFFType& forceFields) {

    //init all forcefields
    _FFManager._forceFields.push_back(
        new FilamentFF(forceFields.FStretchingType,
                       forceFields.FBendingType,
                       forceFields.FTwistingType));
    
    _FFManager._forceFields.push_back(
        new LinkerFF(forceFields.LStretchingType,
                     forceFields.LBendingType,
                     forceFields.LTwistingType));
    
    _FFManager._forceFields.push_back(
        new MotorGhostFF(forceFields.MStretchingType,
                         forceFields.MBendingType,
                         forceFields.MTwistingType));
    
    _FFManager._forceFields.push_back(
        new BranchingFF(forceFields.BrStretchingType,
                        forceFields.BrBendingType,
                        forceFields.BrDihedralType,
                        forceFields.BrPositionType));
    
    _FFManager._forceFields.push_back(
        new VolumeCylindricalFF(forceFields.VolumeFFType));
    
    _FFManager._forceFields.push_back(
        new BoundaryFF(forceFields.BoundaryFFType));

}

