
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

#ifndef M3SYM_MController_h
#define M3SYM_MController_h

#include <vector>

#include "common.h"

#include "Minimizer.h"
#include "ForceFieldManager.h"

#include "Parser.h"

//FORWARD DECLARATIONS
class SubSystem;
class PolakRibiere;
class FletcherRieves;

/// Used to initialize, control, and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize [ForceFields]
 *  (@ref ForceField), given an initial selection of which force fields should be 
 *  included. It can compute forces and energies over all force fields, as well as run
 *  energy [Minimization] (@ref Minimizer) algorithms.
 */
class MController {
    
private:
    ForceFieldManager _FFManager;  ///< Container and methods for all force
                                   ///< fields in system
    vector<Minimizer*> _minimizerAlgorithms; ///< Vector with algorithms
                                             ///< for system minimization
    
    SubSystem* _subSystem; ///< A pointer to the subsystem

    /// Initialize the force-fields used in the simualtion
    void initializeFF (MechanicsFFType& forceFields);
    
    /// Initialize the minimization algorithms used in the simulation
    void initializeMinAlgorithms (MechanicsAlgorithm& Minimizers);

public:
    /// Constructor which sets a subsystem pointer
    MController(SubSystem* s) {_subSystem = s;}
    
    /// Initialze the force fields and minimizers used
    void initialize(MechanicsFFType forceFields, MechanicsAlgorithm Minimizers) {
        initializeFF(forceFields);
        initializeMinAlgorithms(Minimizers);
    }

    /// Run a minimization on the system using the chosen algorithm
    void run() {  _minimizerAlgorithms[0]->equlibrate(_FFManager); }
    
};

#endif
