
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

///FORWARD DECLARATIONS
class SubSystem;

///Controller class is used to initialize, manage, and run an entire simulation

/*!
 *
 *
 *
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

    ///Update the system, called in run
    void updateSystem();
    
public:
    ///Default constructor and destructor
    Controller(SubSystem* s) : _mController(s), _cController(s), _gController(), _subSystem(s) { };
    ~Controller() {};
    
    ///Initialize the system, given an input and output directory
    void initialize(string inputDirectory, string outputDirectory);
    ///Run the simulation
    void run();
    
};

#endif
