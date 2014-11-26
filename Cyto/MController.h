
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

#ifndef M3SYM_MController_h
#define M3SYM_MController_h
#include <vector>

#include "common.h"

#include "Minimizer.h"
#include "ForceFieldManager.h"
#include "Parser.h"

///FORWARD DECLARATIONS
class SubSystem;
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
    ForceFieldManager _FFManager;  ///< container and methods for all force fields in system
    vector<Minimizer*> _minimizerAlgorithms; ///<vector with algorythms for system equlibration
    SubSystem* _subSystem; ///< ptr to the subsystem
    
    ///Initialize the MController using a list of vector names
    ///@param forceFields - a list of forcefields to be added
    void initializeFF (MechanicsFFType& forceFields);
    
    
    void initializeMinAlgorithms (MechanicsAlgorithm& Minimizers);
    
    
public:
    MController(SubSystem* s) {_subSystem = s;}
    
    void initialize(MechanicsFFType forceFields, MechanicsAlgorithm Minimizers)
    {
        initializeFF(forceFields);
        initializeMinAlgorithms(Minimizers);
    }

    ///Run minimization on the system using the chosen algorithm
    void run() {
        _minimizerAlgorithms[0]->equlibrate(_FFManager);
    }
    
};

#endif
