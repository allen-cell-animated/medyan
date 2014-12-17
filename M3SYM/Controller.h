
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

#include "MController.h"
#include "GController.h"
#include "CController.h"

//FORWARD DECLARATIONS
class SubSystem;

/// Used to initialize, manage, and run an entire simulation.

/*!
 *  The Controller is initialized in the main program, and initializes the SubSystem given an initial input directory.
 *  After initialization of all member controllers, the Controller class can run a simulation given the already read
 *  input parameters, by iteratively switching between mechanical equilibration and stochastic chemical steps.
 */

class Controller {

private:
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls

    MController _mController; ///< Chemical controller used
    CController _cController; ///< Mechanical controller used
    GController _gController; ///< Geometry controller used
    
    string _inputDirectory; ///< Input directory being used
    string _outputDirectory; ///< Output directory being used
    
    int _numSteps; ///< Number of chemical steps we are running
    int _numStepsPerMech; ///< Number of chemical steps before mechanical equilibration
    int _numStepsPerSnapshot; ///< Number of steps before a snapshot is recorded
    int _numStepsPerNeighbor; ///< NUmber of steps before a neighbor list update

    /// Update the system, called in run. Will update all positions and reactions
    void updateSystem();
    
    /// Update neighbors lists, called in run
    void updateNeighborLists();
    
public:
    Controller(SubSystem* s) : _mController(s), _cController(s), _gController(), _subSystem(s) {};
    ~Controller() {};
    
    ///Initialize the system, given an input and output directory
    void initialize(string inputDirectory, string outputDirectory);
    ///Run the simulation
    void run();
    
};

#endif
