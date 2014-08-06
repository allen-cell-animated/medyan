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
#include "MController.h"
#include "GController.h"
#include "CController.h"

class SubSystem;

///Controller is used to initialize, manage, and run an entire simulation

class Controller {
    
private:
    
    SubSystem* _subSystem; /// SubSystem that this controller is in
    
    GController _gController; ///< Geometry Controller
    MController _mController; ///< Chemical Controller
    CController _cController; ///< Mechanical Controller
    
public:
    
    void initialize() {
        
        ///Call all controllers to initialize
        _gController.initialize(1, {10,10,10}, {1000.0,1000.0,1000.0});
        _cController.initialize("NRM", "Simple");
        _mController.initialize({""});
        
        //create filaments here
        
        
    }
    
    void run() {
        
        while(true) {
            _cController.run(1000);
            _mController.ComputeEnergy();
            _mController.ComputeForces();
            _mController.run(_subSystem, "FletcherRieves")
        }
        
    }
    
    
    
};


#endif /* defined(__Cyto__Controller__) */
