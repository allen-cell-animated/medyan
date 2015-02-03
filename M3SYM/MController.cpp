
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

#include "MController.h"

#include "SubSystem.h"

#include "FilamentFF.h"
#include "LinkerFF.h"
#include "MotorGhostFF.h"
#include "BoundaryFF.h"
#include "VolumeCylindricalFF.h"
#include "BranchingFF.h"

#include "ConjugateGradient.h"

void MController::initializeFF (MechanicsFFType& forceFields) {

    cout << endl;
    _FFManager._forceFields.push_back(
        new FilamentFF(forceFields.FStretchingType,
                       forceFields.FBendingType,
                       forceFields.FTwistingType));
    
    cout << "Filament force field initialized: " <<endl;
    if(forceFields.FStretchingType != "")
        cout << "Stretching   : " << forceFields.FStretchingType << endl;
    if(forceFields.FBendingType != "")
        cout << "Bending      : " << forceFields.FBendingType<< endl;
    if(forceFields.FTwistingType != "")
        cout << "Twisting     : " << forceFields.FTwistingType <<endl;
    
    _FFManager._forceFields.push_back(
        new LinkerFF(forceFields.LStretchingType,
                     forceFields.LBendingType,
                     forceFields.LTwistingType));
    
    cout << "Linker force field initialized:"<<endl;
    if(forceFields.LStretchingType != "")
        cout << "Stretching   : " << forceFields.LStretchingType<< endl;
    if(forceFields.LBendingType != "")
        cout << "Bending      : " << forceFields.LBendingType << endl;
    if(forceFields.LTwistingType != "")
        cout << "Twisting     : " << forceFields.LTwistingType <<endl;
    
    _FFManager._forceFields.push_back(
        new MotorGhostFF(forceFields.MStretchingType,
                         forceFields.MBendingType,
                         forceFields.MTwistingType));
    
    cout << "Motor force field initialized:"<<endl;
    if(forceFields.MStretchingType != "")
        cout << "Stretching   : " << forceFields.MStretchingType<< endl;
    if(forceFields.MBendingType != "")
        cout << "Bending      : " << forceFields.MBendingType<< endl;
    if(forceFields.MTwistingType != "")
        cout << "Twisting     : " << forceFields.MTwistingType <<endl;
    
    _FFManager._forceFields.push_back(
        new BranchingFF(forceFields.BrStretchingType,
                        forceFields.BrBendingType,
                        forceFields.BrDihedralType,
                        forceFields.BrPositionType));
    
    cout << "Branching force field initialized:"<<endl;
    if(forceFields.BrStretchingType != "")
        cout << "Stretching   : " << forceFields.BrStretchingType<< endl;
    if(forceFields.BrBendingType != "")
        cout << "Bending      : " << forceFields.BrBendingType<< endl;
    if(forceFields.BrDihedralType != "")
        cout << "Dihedral     : " << forceFields.BrDihedralType <<endl;
    if(forceFields.BrPositionType != "")
        cout << "Position     : " << forceFields.BrPositionType <<endl;
    
    _FFManager._forceFields.push_back(
        new VolumeCylindricalFF(forceFields.VolumeFFType));
    
    cout << "Volume force field initialized: " <<endl;
    if(forceFields.VolumeFFType != "")
        cout << "Volume       : " << forceFields.VolumeFFType << endl;
    
    _FFManager._forceFields.push_back(
        new BoundaryFF(forceFields.BoundaryFFType, "", ""));
    
    cout << "Boundary force field initialized: " <<endl;
    if(forceFields.BoundaryFFType != "")
        cout << "Boundary     : " << forceFields.BoundaryFFType << endl;
    
}

void MController::initializeMinAlgorithms (MechanicsAlgorithm& MAlgorithm) {
    
    if (MAlgorithm.ConjugateGradient == "FLETCHERRIEVES")
        _minimizerAlgorithms.push_back(
            new ConjugateGradient<FletcherRieves>(MAlgorithm.gradientTolerance));
    
    else if (MAlgorithm.ConjugateGradient == "POLAKRIBIERE")
        _minimizerAlgorithms.push_back(
            new ConjugateGradient<PolakRibiere>(MAlgorithm.gradientTolerance));
    
    else {
        cout << "Conjugate gradient method not yet implemented. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
}

