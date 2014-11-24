//
//  Controller.h
//  Cyto
//
//  Created by James Komianos on 8/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Controller__
#define __Cyto__Controller__

#include <iostream>
#include "common.h"

#include "MController.h"
#include "GController.h"
#include "CController.h"

///FORWARD DECLARATIONS
class SubSystem;

///Controller is used to initialize, manage, and run an entire simulation
class Controller {

private:
    SubSystem *_subSystem; ///< SubSystem that this controller is in

    MController _mController; ///< Chemical Controller
    CController _cController; ///< Mechanical Controller
    GController _gController; ///< Geometry controller
    
    string _inputDirectory; ///< input directory being used
    string _outputDirectory; ///< output directory being used
    
    int _numSteps; ///< number of chemical steps we are running
    int _numStepsPerMech; ///<number of chemical steps before mechanical equilibration

public:
    ///Default constructor and destructor
    Controller(SubSystem* s) : _mController(s), _cController(s), _gController(), _subSystem(s) { };
    ~Controller() {};
    
    ///Initialize the system
    void initialize(string inputDirectory, string outputDirectory);
    ///update positions of all elements in system
    void updatePositions();
    ///Run the simulation
    void run();
    
};

#endif /* defined(__Cyto__Controller__) */
