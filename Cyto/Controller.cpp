//
//  Controller.cpp
//  Cyto
//
//  Created by James Komianos on 9/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Controller.h"

void Controller::initialize(std::string inputFile) {
    
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
    
#ifdef CHEMISTRY
    ///Initialize chemical controller
    std::cout << "Initializing chemistry...";
    ///read algorithm
    CAlgorithm = p.readChemistryAlgorithm();
    ChemistrySetup CSetup = p.readChemistrySetup();
    
    //num steps for sim
    _numSteps = CAlgorithm.numSteps;
    _numStepsPerMech = CAlgorithm.numStepsPerMech;
    
    ChemistrySpeciesAndReactions chemSR;
    
    if(CSetup.inputFile != "") {
        ChemistryParser cp(CSetup.inputFile);
        chemSR = cp.readChemistryInput();
    }
    else {
        std::cout << "Need to specify a chemical input file. Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }
    _cController.initialize(CAlgorithm.algorithm, "SIMPLE", chemSR);
    ChemSim::printReactions();
    std::cout << "Done." <<std::endl;
#endif
#ifdef MECHANICS
    ///Initialize Mechanical controller
    std::cout << "Initializing mechanics...";
    _mController.initialize(MTypes, MAlgorithm);
    std::cout << "Done." <<std::endl;
    
    ///Initialize boundary
    std::cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->AddBoundary(new BoundaryCubic());
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
    
    
    //ChemSim::printReactions();
}

void Controller::run() {
    
    std::chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = std::chrono::high_resolution_clock::now();
    ///Set up filament output file
    Output o("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/filamentoutput.txt");
    o.printBasicSnapshot(0);
    
#if defined(CHEMISTRY)
    for(int i = 0; i < _numSteps; i+=_numStepsPerMech) {
        _cController.run(_numStepsPerMech);
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
        _mController.run();
        o.printBasicSnapshot(i + _numStepsPerMech);
#elif defined(MECHANICS)
        _mController.run();
        o.printBasicSnapshot(1);
#else
        o.printBasicSnapshot(i);
#endif
#if defined(CHEMISTRY)
    }
#endif
    chk2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_run(chk2-chk1);
    
    std::cout << "Time elapsed for run: dt=" << elapsed_run.count() << std::endl;
    std::cout << "Done with simulation!" << std::endl;
}



