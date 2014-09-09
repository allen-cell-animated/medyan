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
    SubSystem* _subSystem; ///< SubSystem that this controller is in
    
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
        ChemistryParameters CParams;
        BoundaryParameters BParams;
        GeometryParameters GParams;
        MechanicsFFType MTypes;
        MechanicsParameters MParams;
        
        ///read if activated
        if(_mechanics) {
            MParams = p.readMechanicsParameters();
            MTypes = p.readMechanicsFFType();
            BParams = p.readBoundaryParameters();
            
        }
        if(_chemistry) {
            CParams = p.readChemistryParameters();
        }
        
        ///Always read geometry
        GParams = p.readGeometryParameters();
        
        ///Check input
        if(!p.checkInput(CParams, BParams, GParams, MTypes, MParams)) exit(EXIT_FAILURE);
        
        ///CALLING ALL CONTROLLERS TO INITIALIZE
        
        ///Construct grid and compartment vector, initialize grid
        std::vector<int> grid = {};
        std::vector<double> compartmentSize = {};
    
        ///Since we only have one subsystem, this is simple.
        ///In the future, controller will divide up grid into subsystems.
        if(GParams.nDim >= 1) {
            grid.push_back(GParams.NX);
            compartmentSize.push_back(GParams.compartmentSizeX);
        }
        if(GParams.nDim >= 2) {
            grid.push_back(GParams.NY);
            compartmentSize.push_back(GParams.compartmentSizeY);
        }
        if(GParams.nDim == 3) {
            grid.push_back(GParams.NZ);
            compartmentSize.push_back(GParams.compartmentSizeZ);
        }
        GController::initializeGrid(GParams.nDim, grid, compartmentSize);
        
        ///Initialize chemical controller
        if(_chemistry)
            _cController.initialize(CParams.algorithm, CParams.setup);
        
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
