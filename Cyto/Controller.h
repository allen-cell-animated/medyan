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
#include "Parser.h"

class SubSystem;

///Controller is used to initialize, manage, and run an entire simulation
class Controller {
    
private:
    SubSystem _subSystem; ///< SubSystem that this controller is in
    
    MController _mController; ///< Chemical Controller
    CController _cController; ///< Mechanical Controller
    
    bool _mechanics; ///< are we running mechanics?
    bool _chemistry; ///< are we running chemistry?
    
public:
    void initialize(std::string inputFile) {
        
        ///Parse input, get parameters
        Parser p(inputFile);
        
        _mechanics = p.mechanics();
        _chemistry = p.chemistry();
        
        ///Parameters for input
        ChemistryAlgorithm CAlgorithm;
        MechanicsAlgorithm MAlgorithm;
        MechanicsFFType MTypes;
        MechanicsParameters MParams;
        BoundaryParameters BParams;
        GeometryParameters GParams;
        
        ///read if activated
        if(_mechanics) {
            MAlgorithm = p.readMechanicsAlgorithm();
            MParams = p.readMechanicsParameters();
            MTypes = p.readMechanicsFFType();
            BParams = p.readBoundaryParameters();
            
        }
        if(_chemistry) {
            CAlgorithm = p.readChemistryAlgorithm();
        }
        ///Always read geometry
        GParams = p.readGeometryParameters();
        
        ///Check input
        if(!p.checkInput(CAlgorithm, MAlgorithm, MTypes, MParams, BParams, GParams)) exit(EXIT_FAILURE);
        
        ///CALLING ALL CONTROLLERS TO INITIALIZE
        GController::initializeGrid(GParams.nDim, {GParams.NX, GParams.NY, GParams.NZ},
                                    {GParams.compartmentSizeX, GParams.compartmentSizeY, GParams.compartmentSizeZ});
        
        ///Initialize chemical controller
        if(_chemistry)
            _cController.initialize(CAlgorithm.algorithm, CAlgorithm.setup);
        
        ///Initialize Mechanical controller
        if(_mechanics) {
            
            
            
            _mController.initialize({""});
        }
        else {
            ///initialize mechanical controller with no forcefields?
            
            
        }
        
        //create filaments here
        
        
    }
    
    void run() {
        
        while(true) {
            
        }
        
    }
    
    
    
};


#endif /* defined(__Cyto__Controller__) */
