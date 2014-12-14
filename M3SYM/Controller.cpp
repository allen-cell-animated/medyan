
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

#include <random>

#include "Controller.h"

#include "Parser.h"
#include "Output.h"
#include "SubSystem.h"
#include "BoundaryImpl.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

void Controller::initialize(string inputDirectory, string outputDirectory) {
    
    cout << "******************** M3SYM **********************" << endl;
    
    cout.precision(10);
    
    //init input directory
    _inputDirectory = inputDirectory;
    _outputDirectory = outputDirectory;
    
    string inputFile = inputDirectory + "testsysteminput.txt";
    
    //Parse input, get parameters
    SystemParser p(inputFile);
    
#ifdef MECHANICS
    //read algorithm and types
    auto MTypes = p.readMechanicsFFType(); auto MAlgorithm = p.readMechanicsAlgorithm();
    
    //read const parameters
    p.readMechanicsParameters();
#endif
    //Always read boundary type
    auto BTypes = p.readBoundaryType(); p.readBoundaryParameters();
    
    //Always read geometry
    p.readGeometryParameters();
    
    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "Initializing geometry...";
    _gController.initializeGrid();
    cout << "Done." << endl;
    
#ifdef MECHANICS
    //Initialize Mechanical controller
    cout << "Initializing mechanics...";
    _mController.initialize(MTypes, MAlgorithm);
    cout << "Done." <<endl;
    
#endif
    //Initialize boundary
    cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->addBoundary(new BoundaryCubic());
    }
    else if(BTypes.boundaryShape == "SPHERICAL") {
        _subSystem->addBoundary(new BoundarySpherical(SystemParameters::Boundaries().diameter));
    }
    else if(BTypes.boundaryShape == "CAPSULE") {
        _subSystem->addBoundary(new BoundaryCapsule(SystemParameters::Boundaries().diameter));
    }
    else{
        cout << endl << "Given boundary not yet implemented. Exiting" <<endl;
        exit(EXIT_FAILURE);
    }
    cout << "Done." <<endl;
    
#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController.activateCompartments(_subSystem->getBoundary());
    
    //read parameters
    p.readChemistryParameters();
    
    //Initialize chemical controller
    cout << "Initializing chemistry...";
    //read algorithm
    auto CAlgorithm = p.readChemistryAlgorithm(); auto CSetup = p.readChemistrySetup();
    
    //num steps for sim
    _numSteps = CAlgorithm.numSteps;
    _numStepsPerMech = CAlgorithm.numStepsPerMech;
    _numStepsPerSnapshot = CAlgorithm.numStepsPerSnapshot;
    
    ChemistryData chem;
    
    if(CSetup.inputFile != "") {
        ChemistryParser cp(_inputDirectory + CSetup.inputFile);
        chem = cp.readChemistryInput();
    }
    else {
        cout << "Need to specify a chemical input file. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    _cController.initialize(CAlgorithm.algorithm, "", chem);
    cout << "Done." << endl;
#endif
    
    //Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
    vector<vector<vector<double>>> filamentData;
    cout << "Initializing filaments...";
    
    if(FSetup.inputFile != "") {
        FilamentParser fp(_inputDirectory + FSetup.inputFile);
        filamentData = fp.readFilaments();
    }
    else {
        //Create random distribution of filaments
        default_random_engine generator;
        uniform_real_distribution<double> dU(0.0,1.0);
        uniform_real_distribution<double> dUNeg(-1,1);
        
        int filamentCounter = 0;
        while (filamentCounter < FSetup.numFilaments) {
            
            double firstX = dU(generator) * SystemParameters::Geometry().compartmentSizeX *
                                            SystemParameters::Geometry().NX;
            double firstY = dU(generator) * SystemParameters::Geometry().compartmentSizeY *
                                            SystemParameters::Geometry().NY;
            double firstZ = dU(generator) * SystemParameters::Geometry().compartmentSizeZ *
                                            SystemParameters::Geometry().NZ;
            
            double directionX = dUNeg(generator);
            double directionY = dUNeg(generator);
            double directionZ = dUNeg(generator);
            
            //Create a random filament vector one cylinder long
            vector<double> firstPoint = {firstX, firstY, firstZ};
            
            auto normFactor = sqrt(directionX * directionX +
                                   directionY * directionY +
                                   directionZ * directionZ);
            
            vector<double> direction = {directionX/normFactor,
                                        directionY/normFactor,
                                        directionZ/normFactor};
            
            vector<double> secondPoint = nextPointProjection(firstPoint,
               (double)SystemParameters::Geometry().cylinderSize - 0.01, direction);
            
            if(_subSystem->getBoundary()->within(firstPoint)
               && _subSystem->getBoundary()->within(secondPoint)) {
                filamentData.push_back({firstPoint, secondPoint});
                filamentCounter++;
            }
        }
    }
    //add filaments
    _subSystem->addNewFilaments(filamentData);
    cout << "Done. " << filamentData.size() << " filaments created." << endl;
    
    //First update of system
    updateSystem();
}

void Controller::updateSystem() {
    
    /// update all reactables and moveables
    for(auto &f : *FilamentDB::instance()) {
       
        for (auto cylinder : f->getCylinderVector()){
            cylinder->getFirstBead()->updatePosition();
            cylinder->updatePosition();
            cylinder->updateReactionRates();
        }
        //update last bead
        f->getCylinderVector().back()->getSecondBead()->updatePosition();
    }
    for(auto &l : *LinkerDB::instance()) {
        l->updatePosition();
        l->updateReactionRates();
    }
    for(auto &m : *MotorGhostDB::instance()) {
        m->updatePosition();
        m->updateReactionRates();
    }
    for(auto &b : *BranchingPointDB::instance()) {
        b->updatePosition();
        b->updateReactionRates();
    }
    
    //reset neighbor lists
    NeighborListDB::instance()->resetAll();
    
#ifdef CHEMISTRY
    //Update cross cylinder reactions
    for(auto &filament : *FilamentDB::instance())
        for (auto &cylinder : filament->getCylinderVector())
            ChemManager::updateCCylinder(cylinder->getCCylinder());
#endif
    
}

void Controller::run() {
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    ///Set up filament output file
    Output o(_outputDirectory + "filamentoutput.txt");
    o.printBasicSnapshot(0);
    
    cout << "Starting simulation..." << endl;
    
    //perform first minimization
#ifdef MECHANICS
    _mController.run();
    updateSystem();
#endif
    
#if defined(CHEMISTRY)
    for(int i = 0; i < _numSteps; i+=_numStepsPerMech) {
        cout << "Current simulation time = "<< tau() << endl;
        //run ccontroller
        if(!_cController.run(_numStepsPerMech)) break;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
        //run mcontroller, update system
        _mController.run();
        updateSystem();
        
        if(i % _numStepsPerSnapshot == 0)
            o.printBasicSnapshot(i + _numStepsPerMech);
#else
        updateSystem();
        if(i % _numStepsPerSnapshot == 0)
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



