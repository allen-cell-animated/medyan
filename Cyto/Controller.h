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
#include "Output.h"

class SubSystem;

///Controller is used to initialize, manage, and run an entire simulation
class Controller {

private:
    SubSystem *_subSystem; ///< SubSystem that this controller is in

    MController _mController; ///< Chemical Controller
    CController _cController; ///< Mechanical Controller
    
    int _numSteps; ///< number of chemical steps we are running
    int _numStepsPerMech; ///<number of chemical steps before mechanical equilibration

public:
    ///Default constructor and destructor
    Controller(SubSystem* s) : _mController(s), _cController(s), _subSystem(s) { };

    ~Controller() {};

    void initialize(std::string inputFile) {
    
        ///Parse input, get parameters
        SystemParser p(inputFile);
        
        ///Parameters for input
        ChemistryAlgorithm CAlgorithm; MechanicsAlgorithm MAlgorithm;
        MechanicsFFType MTypes; BoundaryType BTypes;
        
#ifdef MECHANICS
        ///read algorithm and types
        MTypes = p.readMechanicsFFType(); BTypes = p.readBoundaryType();
        MAlgorithm = p.readMechanicsAlgorithm();

        ///read const parameters
        p.readMechanicsParameters(); p.readBoundaryParameters();
#endif
        ///Always read geometry
        p.readGeometryParameters();

        ///CALLING ALL CONTROLLERS TO INITIALIZE
        ///Initialize geometry controller
        std::cout << "Initializing geometry...";
        GController::initializeGrid();
        std::cout << "Done." << std::endl;

        ///Initialize chemical controller
#ifdef CHEMISTRY
        std::cout << "Initializing chemistry...";
        
        ///read algorithm
        CAlgorithm = p.readChemistryAlgorithm();
        _numSteps = CAlgorithm.numSteps;
        _numStepsPerMech = CAlgorithm.numStepsPerMech;
        
        
        _cController.initialize(CAlgorithm.algorithm, CAlgorithm.setup);
        ChemSim::printReactions();
        std::cout << "Done." <<std::endl;
#endif
#ifdef MECHANICS
        ///Initialize Mechanical controller
        std::cout << "Initializing mechanics...";
        ///Initialize mcontroller
        _mController.initialize(MTypes, MAlgorithm);
        std::cout << "Done." <<std::endl;

        std::cout << "Initializing boundary...";
        ///Initialize boundary
        if(BTypes.boundaryShape == "CUBIC") {
            //_subSystem->AddBoundary(new BoundaryCubic());
        }
        else{
            std::cout << std::endl << "Given boundary not yet implemented. Exiting" <<std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Done." <<std::endl;
#endif

        ///Read filament setup, parse filament input file if needed
        FilamentSetup FSetup = p.readFilamentSetup();
        std::vector<std::vector<std::vector<double>>> filamentData;
        std::cout << "Initializing filaments...";
        
        if(FSetup.inputFile != "") {
            FilamentParser fp(FSetup.inputFile);
            filamentData = fp.readFilaments();
        }
        else {
            std::cout<< std::endl << "Random filament distributions not yet implemented. Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        ///Create filaments
        _subSystem->AddNewFilaments(filamentData);
        std::cout << "Done." <<std::endl;
    }

    void run() {
        
        ///Set up filament output file
        Output o("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/filamentoutput.txt");

#if defined(MECHANICS) && defined(CHEMISTRY)
        for(int i = 0; i < _numSteps; i+=_numStepsPerMech) {
            _cController.run(_numStepsPerMech);
#elif defined(CHEMISTRY)
        _cController.run(_numSteps);
#endif
#if defined(MECHANICS)
        _mController.run();
#endif
        ///Filament output
        o.printFilaments();
#if defined(MECHANICS) && defined(CHEMISTRY)
        }
#endif
        std::cout << "Done with simulation!" << std::endl;
    }
    
};

#endif /* defined(__Cyto__Controller__) */