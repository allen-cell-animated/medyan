
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

#ifndef M3SYM_Controller_h
#define M3SYM_Controller_h

#include "common.h"

#include "Output.h"
#include "MController.h"
#include "GController.h"
#include "CController.h"

//FORWARD DECLARATIONS
class SubSystem;

/// Used to initialize, manage, and run an entire simulation.

/*!
 *  The Controller is initialized in the main program, and initializes the SubSystem 
 *  given an initial input directory. After initialization of all member controllers,
 *  the Controller class can run a simulation given the already read input parameters, 
 *  by iteratively switching between mechanical equilibration and stochastic chemical 
 *  steps.
 */
class Controller {

private:
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls

    MController* _mController; ///< Chemical controller used
    CController* _cController; ///< Mechanical controller used
    GController* _gController; ///< Geometry controller used
    
    string _inputDirectory; ///< Input directory being used
    string _outputDirectory; ///< Output directory being used
    
    vector<Output*> _outputs; ///< Vector of specified outputs
    
    int _numTotalSteps; ///< Number of chemical steps we are running, or specify the time
    double _runTime; ///< Total desired runtime for simulation
    
    int _numStepsPerSnapshot; ///< Number of steps before a snapshot is recorded
    double _snapshotTime; ///< Desired time for each snapshot
    
    int _numChemSteps; ///< Number of consecutive chemical steps
    int _numStepsPerNeighbor; ///< NUmber of steps before a neighbor list update
    
    /// Update the positions of all elements in the system
    void updatePositions();
    
#ifdef DYNAMICRATES
    /// Update the reaction rates of all elements in the system
    void updateReactionRates();
#endif
    
    /// Update neighbors lists, called in run
    void updateNeighborLists(bool updateReactions=false);
    
public:
    Controller(SubSystem* s);
    ~Controller() {};
    
    ///Initialize the system, given an input and output directory
    void initialize(string inputDirectory, string outputDirectory);
    ///Run the simulation
    void run();
    
};

#endif
