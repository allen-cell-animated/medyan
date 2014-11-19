//
//  Controller.cpp
//  Cyto
//
//  Created by James Komianos on 9/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Controller.h"

#include "Parser.h"
#include "Output.h"

#include "BoundaryImpl.h"
#include "NeighborListDB.h"

#include "MathFunctions.h"

using namespace mathfunc;

void Controller::initialize(string inputFile) {
    
    cout << "******************** CYTOSIM **********************" << endl;
    
    ///Parse input, get parameters
    SystemParser p(inputFile);
    
    ///Parameters for input
    ChemistryAlgorithm CAlgorithm; MechanicsAlgorithm MAlgorithm;
    MechanicsFFType MTypes; BoundaryType BTypes;
    
#ifdef MECHANICS
    ///read algorithm and types
    MTypes = p.readMechanicsFFType(); MAlgorithm = p.readMechanicsAlgorithm();
    
    ///read const parameters
    p.readMechanicsParameters();
#endif
    ///Always read boundary type
    BTypes = p.readBoundaryType(); p.readBoundaryParameters();
    
    ///Always read geometry
    p.readGeometryParameters();
    
    ///CALLING ALL CONTROLLERS TO INITIALIZE
    ///Initialize geometry controller
    cout << "Initializing geometry...";
    _gController.initializeGrid();
    cout << "Done." << endl;
    
    ///Initialize boundary
    cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->addBoundary(new BoundaryCubic());
    }
    else if(BTypes.boundaryShape == "SPHERICAL") {
        _subSystem->addBoundary(new BoundarySpherical());
    }
    else{
        cout << endl << "Given boundary not yet implemented. Exiting" <<endl;
        exit(EXIT_FAILURE);
    }
    cout << "Done." <<endl;
    
#ifdef CHEMISTRY
    ///Activate necessary compartments for diffusion
    _gController.activateCompartments(_subSystem->getBoundary());
    
    ///read parameters
    p.readChemistryParameters();
    
    ///Initialize chemical controller
    cout << "Initializing chemistry...";
    ///read algorithm
    CAlgorithm = p.readChemistryAlgorithm(); ChemistrySetup CSetup = p.readChemistrySetup();
    
    //num steps for sim
    _numSteps = CAlgorithm.numSteps;
    _numStepsPerMech = CAlgorithm.numStepsPerMech;
    
    ChemistryData chem;
    
    if(CSetup.inputFile != "") {
        ChemistryParser cp(CSetup.inputFile);
        chem = cp.readChemistryInput();
    }
    else {
        cout << "Need to specify a chemical input file. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    _cController.initialize(CAlgorithm.algorithm, "", chem);
    cout << "Done." <<endl;
#endif
#ifdef MECHANICS
    ///Initialize Mechanical controller
    cout << "Initializing mechanics...";
    _mController.initialize(MTypes, MAlgorithm);
    cout << "Done." <<endl;
    
#endif
    
    ///Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
    vector<vector<vector<double>>> filamentData;
    cout << "Initializing filaments...";
    
    if(FSetup.inputFile != "") {
        FilamentParser fp(FSetup.inputFile);
        filamentData = fp.readFilaments();
    }
    else {
        ///Create random distribution of filaments
        default_random_engine generator;
        uniform_real_distribution<double> dU(0,1);
        uniform_real_distribution<double> dUNeg(-1,1);
        
        for(int i = 0; i < FSetup.numFilaments; i++) {
            
            double firstX = dU(generator) * SystemParameters::Geometry().compartmentSizeX *
                                            SystemParameters::Geometry().NX;
            double firstY = dU(generator) * SystemParameters::Geometry().compartmentSizeY *
                                            SystemParameters::Geometry().NY;
            double firstZ = dU(generator) * SystemParameters::Geometry().compartmentSizeZ *
                                            SystemParameters::Geometry().NZ;
            
            double directionX = dUNeg(generator);
            double directionY = dUNeg(generator);
            double directionZ = dUNeg(generator);
            
            ///Create a random filament vector one cylinder long
            vector<double> firstPoint = {firstX, firstY, firstZ};
            auto normFactor = sqrt(directionX * directionX +
                                   directionY * directionY +
                                   directionZ * directionZ);
            
            vector<double> direction = {directionX/normFactor,
                                        directionY/normFactor,
                                        directionZ/normFactor};
            vector<double> secondPoint = NextPointProjection(firstPoint,
                    (double)SystemParameters::Geometry().cylinderSize - 0.01, direction);
            
            if(_subSystem->getBoundary()->within(firstPoint)
               && _subSystem->getBoundary()->within(secondPoint)) {
                filamentData.push_back({firstPoint, secondPoint});
            }
        }
    }
    ///add filaments
    _subSystem->addNewFilaments(filamentData);
    cout << "Done. " << filamentData.size() << " filaments created." << endl;
    
    ///First update of positions
    updatePositions();
}

void Controller::updatePositions() {
    
    ///Update bead-boundary interactions
    for(auto b : *BeadDB::instance(BeadDBKey())) b->updatePosition();
    ///Update cylinder positions
    for(auto &c : *CylinderDB::instance(CylinderDBKey())) c->updatePosition();
    ///Update linker positions
    for(auto &l : *LinkerDB::instance(LinkerDBKey())) l->updatePosition();
    ///update motor positions
    for(auto &m : *MotorGhostDB::instance(MotorGhostDBKey())) m->updatePosition();
    
    ///reset neighbor lists
    NeighborListDB::instance(NeighborListDBKey())->resetAll();
    
#ifdef CHEMISTRY
    //Update filament reactions
    for(auto &c : *CylinderDB::instance(CylinderDBKey()))
        ChemManager::updateCCylinder(ChemManagerCylinderKey(), c->getCCylinder());
#endif
    
}

void Controller::run() {
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    ///Set up filament output file
    //Output o("/Users/Konstantin/Documents/Codes/Cyto/CytoRepo/Cyto/filamentoutput.txt");
    Output o("/Users/jameskomianos/Code/CytoSim-Repo/Cyto/filamentoutput.txt");
    o.printBasicSnapshot(0);
    
    cout << "Starting simulation..." << endl;
    
#if defined(CHEMISTRY)
    for(int i = 0; i < _numSteps; i+=_numStepsPerMech) {
        
        cout << "Current simulation time = "<< tau() << endl;
        _cController.run(_numStepsPerMech);
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
        _mController.run();
        updatePositions();
        o.printBasicSnapshot(i + _numStepsPerMech);
#elif defined(MECHANICS)
        _mController.run();
        updatePositions();
        o.printBasicSnapshot(1);
#else
        updatePositions();
        o.printBasicSnapshot(i + _numStepsPerMech);
#endif
        
#if defined(CHEMISTRY)
    }
#endif
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
}



