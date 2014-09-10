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
#include "BoundaryImpl.h"

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
        SystemParser p(inputFile);
        
        _mechanics = p.mechanics();
        _chemistry = p.chemistry();
        
        ///Parameters for input
        ChemistryAlgorithm CAlgorithm; MechanicsAlgorithm MAlgorithm;
        MechanicsFFType MTypes; BoundaryType BTypes;
        
        ///read if activated
        if(_mechanics) {
            ///read algorithm and types
            MTypes = p.readMechanicsFFType();
            BTypes = p.readBoundaryType();
            MAlgorithm = p.readMechanicsAlgorithm();
            
            ///read const parameters
            p.readMechanicsParameters();
            p.readBoundaryParameters();
        }
        if(_chemistry) {
            ///read algorithm
            CAlgorithm = p.readChemistryAlgorithm();
        }
        ///Always read geometry
        p.readGeometryParameters();

        ///CALLING ALL CONTROLLERS TO INITIALIZE
        ///Initialize geometry controller
        short nDim = SystemParameters::Geometry().nDim;
        std::vector<int> grid = {SystemParameters::Geometry().NX,
            SystemParameters::Geometry().NY,
            SystemParameters::Geometry().NZ};
        std::vector<double> compartmentSize = {SystemParameters::Geometry().compartmentSizeX,
            SystemParameters::Geometry().compartmentSizeY,
            SystemParameters::Geometry().compartmentSizeZ};

        GController::initializeGrid(nDim, grid, compartmentSize);
        
        ///Initialize chemical controller
        if(_chemistry)
            _cController.initialize(CAlgorithm.algorithm, CAlgorithm.setup);
        
        ///Initialize Mechanical controller
        if(_mechanics) {
            
            ///Boundary
            if(BTypes.boundaryShape == "CUBIC")
                _subSystem.AddBoundary(new BoundaryCubic());
            else {
                std::cout << "Boundary not yet implemented. Exiting." <<std::endl;
                exit(EXIT_FAILURE);
            }
            
            
            
            _mController.initialize({""});
        }
        else {
            ///initialize mechanical controller with no forcefields?
            
            
        }
        
        ///Read filament setup, parse filament input file if needed
        FilamentSetup FSetup = p.readFilamentSetup();
        std::vector<std::vector<std::vector<double>>> filamentData;
        
        if(FSetup.inputFile != "") {
            FilamentParser fp(FSetup.inputFile);
            filamentData = fp.readFilaments();
        }
        else {
            std::cout << "Random filament distributions not yet implemented. Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        ///Create 
        
    }
    
    void run() {
        
        while(true) {
            
        }
        
    }
    
    
    
};


#endif /* defined(__Cyto__Controller__) */
