
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MController.h"

#include "SubSystem.h"

#include "FilamentFF.h"
#include "LinkerFF.h"
#include "MotorGhostFF.h"
#include "BoundaryFF.h"
#include "BranchingFF.h"
#include "BubbleFF.h"
#include "CylinderVolumeFF.h"

#include "ConjugateGradient.h"

void MController::initializeMinAlgorithms (MechanicsAlgorithm& MAlgorithm) {

    if (MAlgorithm.ConjugateGradient == "FLETCHERRIEVES")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<FletcherRieves>(MAlgorithm.gradientTolerance,
                                              MAlgorithm.maxDistance,
                                              MAlgorithm.lambdaMax,
                                              MAlgorithm.lambdarunningaverageprobability));
    else if (MAlgorithm.ConjugateGradient == "POLAKRIBIERE")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<PolakRibiere>(MAlgorithm.gradientTolerance,
                                            MAlgorithm.maxDistance,
                                            MAlgorithm.lambdaMax,
                                            MAlgorithm.lambdarunningaverageprobability));
    else if (MAlgorithm.ConjugateGradient == "STEEPESTDESCENT")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<SteepestDescent>(MAlgorithm.gradientTolerance,
                                               MAlgorithm.maxDistance,
                                               MAlgorithm.lambdaMax,
                                               MAlgorithm.lambdarunningaverageprobability));

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

    //These FF's have a neighbor list associated with them
    //add to the subsystem's database of neighbor lists.
    auto volumeFF = new CylinderVolumeFF(forceFields.VolumeFFType);
    _FFManager._forceFields.push_back(volumeFF);

    //Get the force field access to the HNLID
    for(auto nl : volumeFF->getNeighborLists()) {

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        //Original neighborlist is created even when HybridNeighborList is use as the
        // neighborLists are implemented as predominantly modifications in the back-end.
        // Front end functions remain the same and function calls are internally
        // redirected to HybridNeighborList functions.
            volumeFF->setHNeighborLists(_subSystem->getHNeighborList());
#endif
        if(nl != nullptr)
            _subSystem->addNeighborList(nl);
    }

    auto boundaryFF = new BoundaryFF(forceFields.BoundaryFFType);
    _FFManager._forceFields.push_back(boundaryFF);
    for(auto nl : boundaryFF->getNeighborLists()) {

        if(nl != nullptr)
#if defined(NLSTENCILLIST) || defined(NLORIGINAL)
            _subSystem->addNeighborList(nl);
#endif
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        _subSystem->addBNeighborList(nl);
#endif
    }

    auto bubbleFF = new BubbleFF(forceFields.BubbleFFType,
                                 forceFields.MTOCFFType);
    _FFManager._forceFields.push_back(bubbleFF);
    for(auto nl : bubbleFF->getNeighborLists()) {

        if(nl != nullptr)
            _subSystem->addBNeighborList(nl);
    }
}
