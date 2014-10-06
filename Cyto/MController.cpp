//
//  MController.cpp
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MController.h"

void MController::initializeFF (MechanicsFFType& forceFields) {
    /// Check if exist!!!
    _FFManager._forceFields.push_back(new FilamentFF(forceFields.FStretchingType, forceFields.FBendingType, forceFields.FTwistingType) );
    std::cout << "Filament force field initialized: " <<std::endl;
    if(forceFields.FStretchingType != "") std::cout << "Stretching: " << forceFields.FStretchingType << std::endl;
    if(forceFields.FBendingType != "") std::cout << "Bending: " << forceFields.FBendingType<< std::endl;
    if(forceFields.FTwistingType != "")std::cout << "Twisting: " << forceFields.FTwistingType <<std::endl;
    
    _FFManager._forceFields.push_back(new LinkerFF(forceFields.LStretchingType, forceFields.LBendingType, forceFields.LTwistingType) );
    std::cout << "Linker force field initialized:"<<std::endl;
    if(forceFields.LStretchingType != "") std::cout << "Stretching: " << forceFields.LStretchingType<< std::endl;
    if(forceFields.LBendingType != "") std::cout << "Bending: " << forceFields.LBendingType<< std::endl;
    if(forceFields.LTwistingType != "") std::cout << "Twisting: " << forceFields.LTwistingType <<std::endl;
    
    _FFManager._forceFields.push_back(new MotorGhostFF(forceFields.MStretchingType, forceFields.MBendingType, forceFields.MTwistingType) );
    std::cout << "Motor force field initialized:"<<std::endl;
    if(forceFields.MStretchingType != "") std::cout << "Stretching: " << forceFields.MStretchingType<< std::endl;
    if(forceFields.MBendingType != "") std::cout << "Bending: " << forceFields.MBendingType<< std::endl;
    if(forceFields.MTwistingType != "") std::cout << "Twisting: " << forceFields.MTwistingType <<std::endl;
    
    _FFManager._forceFields.push_back(new BoundaryFF(forceFields.BoundaryFFType, "", "") );
    std::cout << "Boundary force field initialized: " <<std::endl;
    if(forceFields.BoundaryFFType != "") std::cout << "Boundary: " << forceFields.BoundaryFFType << std::endl;
    
}

void MController::initializeMinAlgorithms (MechanicsAlgorithm& Minimizers) {
    if (Minimizers.ConjugateGradient == "FLETCHERRIEVES") {_minimizerAlgorithms.push_back(new ConjugateGradient<FletcherRieves>() );}
    else if (Minimizers.ConjugateGradient == "POLAKRIBIERE") {_minimizerAlgorithms.push_back(new ConjugateGradient<PolakRibiere>() );}
    else {
        std::cout << "Conjugate gradient method not yet implemented. Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void MController::updatePositions() {

    ///Update bead-boundary interactions (VERY INEFFICIENT)
    for(auto b : *BeadDB::Instance(BeadDBKey())) b->updatePosition();
    
    ///Update cylinder positions (ALSO VERY INEFFICIENT)
    for(auto &c : *CylinderDB::Instance(CylinderDBKey())) c->updatePosition();
    
    ///Update linker positions
    for(auto &l : *LinkerDB::Instance(LinkerDBKey())) l->updatePosition();
    
}

