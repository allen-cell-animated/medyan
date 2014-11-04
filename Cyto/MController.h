//
//  MController.h
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MController__
#define __Cyto__MController__

#include "common.h"
#include "SubSystem.h"
#include "Parser.h"

#include "ForceFieldManager.h"
#include "Minimizer.h"

#include "ForceField.h"
#include "FilamentFF.h"
#include "LinkerFF.h"
#include "VolumeCylindricalFF.h"
#include "BoundaryFF.h"
#include "MotorGhostFF.h"

#include "ConjugateGradient.h"

#include "FilamentDB.h"
#include "CylinderDB.h"
#include "LinkerDB.h"
#include "MotorGhostDB.h"

#include <iostream>
#include <vector>

class PolakRibiere;
class FletcherRieves;


/// MController class is used to initialize and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize force fields, given an initial
 *  selection of which force fields should be included. It can compute forces and energies
 *  over all force fields, as well as run energy minimization algorithms.
 */

class MController {
    
private:
    ForceFieldManager _FFManager;  ///<container and methods for all force fields in system
    std::vector<Minimizer*> _minimizerAlgorithms; ///<vector with algorythms for system equlibration
    SubSystem* _subSystem;
    
    ///Initialize the MController using a list of vector names
    ///@param forceFields - a list of forcefields to be added
    void initializeFF (MechanicsFFType& forceFields);
    void initializeMinAlgorithms (MechanicsAlgorithm& Minimizers);
    
    
public:
    MController(SubSystem* s) {_subSystem = s;}
    
    void initialize(MechanicsFFType forceFields, MechanicsAlgorithm Minimizers )
    {
        initializeFF(forceFields);
        initializeMinAlgorithms(Minimizers);
    }

    ///Run minimization on the system using the chosen algorithm
    void run() {
        _minimizerAlgorithms[0]->Equlibrate(_FFManager);
    }
    
};

#endif /* defined(__Cyto__MController__) */
