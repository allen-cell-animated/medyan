
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
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
#include "CaMKIIingFF.h"
#include "BubbleFF.h"
#include "CylinderVolumeFF.h"

#include "ConjugateGradient.h"

void MController::initializeMinAlgorithms (MechanicsAlgorithm& MAlgorithm) {
    
    if (MAlgorithm.ConjugateGradient == "FLETCHERRIEVES")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<FletcherRieves>(MAlgorithm.gradientTolerance,
                                              MAlgorithm.maxDistance,
                                              MAlgorithm.lambdaMax));
    else if (MAlgorithm.ConjugateGradient == "POLAKRIBIERE")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<PolakRibiere>(MAlgorithm.gradientTolerance,
                                            MAlgorithm.maxDistance,
                                            MAlgorithm.lambdaMax));
    else if (MAlgorithm.ConjugateGradient == "STEEPESTDESCENT")
        _minimizerAlgorithms.push_back(
        new ConjugateGradient<SteepestDescent>(MAlgorithm.gradientTolerance,
                                               MAlgorithm.maxDistance,
                                               MAlgorithm.lambdaMax));
    
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
        new CaMKIIingFF(forceFields.CaMKIIStretchingType,
                        forceFields.CaMKIIBendingType,
                        forceFields.CaMKIIDihedralType,
                        forceFields.CaMKIIPositionType));
    
    //These FF's have a neighbor list associated with them
    //add to the subsystem's database of neighbor lists.
    auto volumeFF = new CylinderVolumeFF(forceFields.VolumeFFType);
    _FFManager._forceFields.push_back(volumeFF);
    for(auto nl : volumeFF->getNeighborLists()) {
        
        if(nl != nullptr)
            _subSystem->addNeighborList(nl);
    }
    
    auto boundaryFF = new BoundaryFF(forceFields.BoundaryFFType);
    _FFManager._forceFields.push_back(boundaryFF);
    for(auto nl : boundaryFF->getNeighborLists()) {
        
        if(nl != nullptr)
            _subSystem->addNeighborList(nl);
    }
    
    auto bubbleFF = new BubbleFF(forceFields.BubbleFFType,
                                 forceFields.MTOCFFType);
    _FFManager._forceFields.push_back(bubbleFF);
    for(auto nl : bubbleFF->getNeighborLists()) {
        
        if(nl != nullptr)
            _subSystem->addNeighborList(nl);
    }
}

