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
    cout << endl;
    _FFManager._forceFields.push_back(new FilamentFF(forceFields.FStretchingType, forceFields.FBendingType, forceFields.FTwistingType) );
    cout << "Filament force field initialized: " <<endl;
    if(forceFields.FStretchingType != "") cout << "Stretching: " << forceFields.FStretchingType << endl;
    if(forceFields.FBendingType != "") cout << "Bending: " << forceFields.FBendingType<< endl;
    if(forceFields.FTwistingType != "")cout << "Twisting: " << forceFields.FTwistingType <<endl;
    
    _FFManager._forceFields.push_back(new LinkerFF(forceFields.LStretchingType, forceFields.LBendingType, forceFields.LTwistingType) );
    cout << "Linker force field initialized:"<<endl;
    if(forceFields.LStretchingType != "") cout << "Stretching: " << forceFields.LStretchingType<< endl;
    if(forceFields.LBendingType != "") cout << "Bending: " << forceFields.LBendingType<< endl;
    if(forceFields.LTwistingType != "") cout << "Twisting: " << forceFields.LTwistingType <<endl;
    
    _FFManager._forceFields.push_back(new MotorGhostFF(forceFields.MStretchingType, forceFields.MBendingType, forceFields.MTwistingType) );
    cout << "Motor force field initialized:"<<endl;
    if(forceFields.MStretchingType != "") cout << "Stretching: " << forceFields.MStretchingType<< endl;
    if(forceFields.MBendingType != "") cout << "Bending: " << forceFields.MBendingType<< endl;
    if(forceFields.MTwistingType != "") cout << "Twisting: " << forceFields.MTwistingType <<endl;
    
    _FFManager._forceFields.push_back(new VolumeCylindricalFF(forceFields.VolumeFFType) );
    cout << "Volume force field initialized: " <<endl;
    if(forceFields.VolumeFFType != "") cout << "Volume: " << forceFields.VolumeFFType << endl;
    
    _FFManager._forceFields.push_back(new BoundaryFF(forceFields.BoundaryFFType, "", "") );
    cout << "Boundary force field initialized: " <<endl;
    if(forceFields.BoundaryFFType != "") cout << "Boundary: " << forceFields.BoundaryFFType << endl;
    
}

void MController::initializeMinAlgorithms (MechanicsAlgorithm& Minimizers) {
    if (Minimizers.ConjugateGradient == "FLETCHERRIEVES") {_minimizerAlgorithms.push_back(new ConjugateGradient<FletcherRieves>() );}
    else if (Minimizers.ConjugateGradient == "POLAKRIBIERE") {_minimizerAlgorithms.push_back(new ConjugateGradient<PolakRibiere>() );}
    else {
        cout << "Conjugate gradient method not yet implemented. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
}

