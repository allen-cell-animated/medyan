
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

#include "FilamentInitializer.h"
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

Controller::Controller(SubSystem* s) : _subSystem(s) {
    
    //init subsystem
    _subSystem = new SubSystem();
    
    //init controllers
    _mController = new MController(_subSystem);
    _cController = new CController(_subSystem);
    _gController = new GController();
}

void Controller::initialize(string inputFile,
                            string inputDirectory,
                            string outputDirectory) {
    
    cout << "******************** M3SYM **********************" << endl;
    
    cout.precision(10);
    
    //init input directory
    _inputDirectory = inputDirectory;
    _outputDirectory = outputDirectory;
    
    //Parse input, get parameters
    SystemParser p(inputFile);
    
    //read output types
    OutputTypes oTypes = p.readOutputTypes();
    
    //snapshot type output
    if(oTypes.basicSnapshot)
        _outputs.push_back(new BasicSnapshot(_outputDirectory + "snapshot.traj"));
    if(oTypes.birthTimes)
        _outputs.push_back(new BirthTimes(_outputDirectory + "birthtimes.traj"));
    if(oTypes.forces)
        _outputs.push_back(new Forces(_outputDirectory + "forces.traj"));
    if(oTypes.stresses)
        _outputs.push_back(new Stresses(_outputDirectory + "stresses.traj"));
    
    
#ifdef MECHANICS
    //read algorithm and types
    auto MTypes = p.readMechanicsFFType();
    auto MAlgorithm = p.readMechanicsAlgorithm();
    
    //read const parameters
    p.readMechanicsParameters();
#endif
    //Always read boundary type
    auto BTypes = p.readBoundaryType();
    p.readBoundaryParameters();
    
    //Always read geometry
    p.readGeometryParameters();
    
    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "Initializing geometry...";
    _gController->initializeGrid();
    cout << "Done." << endl;
    
#ifdef MECHANICS
    //Initialize Mechanical controller
    cout << "Initializing mechanics...";
    _mController->initialize(MTypes, MAlgorithm);
    cout << "Done." <<endl;
    
#endif
    //Initialize boundary
    cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->addBoundary(new BoundaryCubic());
    }
    else if(BTypes.boundaryShape == "SPHERICAL") {
        _subSystem->addBoundary(
            new BoundarySpherical(SystemParameters::Boundaries().diameter));
    }
    else if(BTypes.boundaryShape == "CAPSULE") {
        _subSystem->addBoundary(
            new BoundaryCapsule(SystemParameters::Boundaries().diameter));
    }
    else{
        cout << endl << "Given boundary not yet implemented. Exiting" <<endl;
        exit(EXIT_FAILURE);
    }
    cout << "Done." <<endl;
    
#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController->activateCompartments(_subSystem->getBoundary());
    
    //read parameters
    p.readChemistryParameters();
    
    //Initialize chemical controller
    cout << "Initializing chemistry...";
    //read algorithm
    auto CAlgorithm = p.readChemistryAlgorithm();
    auto CSetup = p.readChemistrySetup();
    
    //num steps for sim
    _numTotalSteps = CAlgorithm.numTotalSteps;
    _runTime = CAlgorithm.runTime;
    
    //if no snapshot step size set, set this to maxint so we use time
    _numStepsPerSnapshot = CAlgorithm.numStepsPerSnapshot;
    if(_numStepsPerSnapshot == 0)
        _numStepsPerSnapshot = numeric_limits<int>::max();
    
    _snapshotTime = CAlgorithm.snapshotTime;
    
    _numChemSteps= CAlgorithm.numChemSteps;
    _numStepsPerNeighbor = CAlgorithm.numStepsPerNeighbor;
    
    ChemistryData chem;
    
    if(CSetup.inputFile != "") {
        ChemistryParser cp(_inputDirectory + CSetup.inputFile);
        chem = cp.readChemistryInput();
    }
    else {
        cout << "Need to specify a chemical input file. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    _cController->initialize(CAlgorithm.algorithm, "", chem);
    cout << "Done." << endl;
#endif
    
#ifdef DYNAMICRATES
    //read dynamic rate parameters
    p.readDynamicRateParameters();
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
        FilamentInitializer* fInit = new RandomFilamentDist();
        filamentData = fInit->createFilaments(_subSystem->getBoundary(),
                                              FSetup.numFilaments,
                                              FSetup.filamentLength);
        delete fInit;
    }
    //add filaments
    _subSystem->addNewFilaments(filamentData);
    cout << "Done. " << filamentData.size() << " filaments created." << endl;
    
    //First update of system
    updatePositions();
    updateNeighborLists();
}

void Controller::updatePositions() {
    
    /// update all moveables
    for(auto &f : *FilamentDB::instance()) {
        for (auto cylinder : f->getCylinderVector()){
            cylinder->getFirstBead()->updatePosition();
            cylinder->updatePosition();
        }
        //update last bead
        f->getCylinderVector().back()->
        getSecondBead()->updatePosition();
    }
    for(auto &l : *LinkerDB::instance())
        l->updatePosition();
    
    for(auto &m : *MotorGhostDB::instance())
        m->updatePosition();
    
    for(auto &b : *BranchingPointDB::instance())
        b->updatePosition();
}

#ifdef DYNAMICRATES
void Controller::updateReactionRates() {
    /// update all reactables
    for(auto &c : _subSystem->getBoundaryCylinders())
        c->updateReactionRates();
    
    for(auto &l : *LinkerDB::instance())
        l->updateReactionRates();
    
    for(auto &m : *MotorGhostDB::instance())
        m->updateReactionRates();
}
#endif

void Controller::updateNeighborLists(bool updateReactions) {
    
    //Full reset of neighbor lists
    NeighborListDB::instance()->resetAll();
    
    if(updateReactions) {
#ifdef CHEMISTRY
        //Update cross cylinder reactions
        for(auto &filament : *FilamentDB::instance())
            for (auto &cylinder : filament->getCylinderVector())
                ChemManager::updateCCylinder(cylinder->getCCylinder());
#endif
    }
}

void Controller::run() {
    
    double tauLastSnapshot = 0;
    double oldTau = 0;
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    ///Print initial configuration
    for(auto o: _outputs) o->print(0);
    
    cout << "Performing an initial minimization..." << endl;
    
    //perform first minimization
#ifdef MECHANICS
    _mController->run();
    
    updatePositions();
    updateNeighborLists(true);
#ifdef DYNAMICRATES
    updateReactionRates();
#endif
    
#endif
    cout << "Starting simulation..." << endl;
    
    //if runtime was specified, use this
    if(_runTime != 0) {
    
#ifdef CHEMISTRY
        int i = 0;
        while(tau() <= _runTime) {
            //run ccontroller
            if(!_cController->run(_numChemSteps)) break;
            
            //add the last step
            tauLastSnapshot += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            _mController->run();
            updatePositions();
            
            if(i % _numStepsPerSnapshot == 0 ||
               tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i + _numChemSteps);
                tauLastSnapshot = 0.0;
            }
#elif defined(CHEMISTRY)
            if(i % _numStepsPerSnapshot == 0
               tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i + _numChemSteps);
                tauLastSnapshot = 0.0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(1);
#endif
            // update neighbor lists
            if(i % _numStepsPerNeighbor == 0 && i != 0)
                updateNeighborLists(true);
#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            i += _numChemSteps;
            oldTau = tau();
        }
#endif
    }
    
    else {
#ifdef CHEMISTRY
        for(int i = 0; i < _numTotalSteps; i+=_numChemSteps) {
            //run ccontroller
            if(!_cController->run(_numChemSteps)) break;
            
            //add the last step
            tauLastSnapshot += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            _mController->run();
            updatePositions();
            
            if(i % _numStepsPerSnapshot == 0 ||
               tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i + _numChemSteps);
                tauLastSnapshot = 0.0;
            }
#elif defined(CHEMISTRY)
            if(i % _numStepsPerSnapshot == 0 ||
               tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i + _numChemSteps);
                tauLastSnapshot = 0.0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(1);
#endif
            // update neighbor lists
            if(i % _numStepsPerNeighbor == 0 && i != 0)
                updateNeighborLists(true);
#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            oldTau = tau();
        }
#endif
    }

    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
}



