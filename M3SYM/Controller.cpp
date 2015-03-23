
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include <random>
#include <chrono>

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

#include "SysParams.h"

Controller::Controller(SubSystem* s) : _subSystem(s) {
    
    //init subsystem
    _subSystem = new SubSystem();
    
    //init controllers
    _mController   = new MController(_subSystem);
    _cController   = new CController(_subSystem);
    _gController   = new GController();
    _drController  = new DRController();
}

void Controller::initialize(string inputFile,
                            string inputDirectory,
                            string outputDirectory) {
    
    //general check of macros
#if defined(DYNAMICRATES) && (!defined(CHEMISTRY) || !defined(MECHANICS))
    cout << "If dynamic rates is turned on, chemistry and mechanics must be "
         << "defined. Please set these compilation macros and try again. Exiting."
         << endl;
    exit(EXIT_FAILURE);
#endif
    
    //init input directory
    _inputDirectory = inputDirectory;
    _outputDirectory = outputDirectory;
    
    //Parse input, get parameters
    SystemParser p(inputFile);
    
    //read output types
    OutputTypes oTypes = p.readOutputTypes();
    
    //snapshot type output
    cout << endl;
    if(oTypes.basicSnapshot)
        _outputs.push_back(new BasicSnapshot(_outputDirectory + "snapshot.traj"));
    if(oTypes.birthTimes)
        _outputs.push_back(new BirthTimes(_outputDirectory + "birthtimes.traj"));
    if(oTypes.forces)
        _outputs.push_back(new Forces(_outputDirectory + "forces.traj"));
    if(oTypes.stresses)
        _outputs.push_back(new Stresses(_outputDirectory + "stresses.traj"));
    
    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters())
        exit(EXIT_FAILURE);
    
    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "---" << endl;
    cout << "Initializing geometry...";
    _gController->initializeGrid();
    cout << "Done." << endl;
    
#ifdef MECHANICS
    //read algorithm and types
    auto MTypes = p.readMechanicsFFType();
    auto MAlgorithm = p.readMechanicsAlgorithm();
    
    //read const parameters
    p.readMechParams();
    
    //Initialize Mechanical controller
    cout << "---" << endl;
    cout << "Initializing mechanics...";
    _mController->initialize(MTypes, MAlgorithm);
    cout << "Done." <<endl;
    
#endif
    //Always read boundary type
    auto BTypes = p.readBoundaryType();
    p.readBoundParams();
    
    //Initialize boundary
    cout << "---" << endl;
    cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->addBoundary(new BoundaryCubic());
    }
    else if(BTypes.boundaryShape == "SPHERICAL") {
        _subSystem->addBoundary(
            new BoundarySpherical(SysParams::Boundaries().diameter));
    }
    else if(BTypes.boundaryShape == "CAPSULE") {
        _subSystem->addBoundary(
            new BoundaryCapsule(SysParams::Boundaries().diameter));
    }
    else{
        cout << endl << "Given boundary not yet implemented. Exiting." <<endl;
        exit(EXIT_FAILURE);
    }
    cout << "Done." <<endl;
    
#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController->activateCompartments(_subSystem->getBoundary());
    
    //read parameters
    p.readChemParams();
    
    //Initialize chemical controller
    cout << "---" << endl;
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
    
    ChemistryData ChemData;
    
    if(CSetup.inputFile != "") {
        ChemistryParser cp(_inputDirectory + CSetup.inputFile);
        ChemData = cp.readChemistryInput();
    }
    else {
        cout << "Need to specify a chemical input file. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    _cController->initialize(CAlgorithm.algorithm, ChemData);
    cout << "Done." << endl;
#endif
    
#ifdef DYNAMICRATES
    cout << "---" << endl;
    cout << "Initializing dynamic rates...";
    //read dynamic rate parameters
    p.readDyRateParams();
    
    //read dynamic rate types
    DynamicRateTypes DRTypes = p.readDynamicRateTypes();
    
    //init controller
    _drController->initialize(DRTypes);
    cout << "Done." << endl;
    
#endif

    //Check consistency of all chemistry and mechanics parameters
    cout << "---" << endl;
    cout << "Checking cross-parameter consistency..." << endl;
#ifdef CHEMISTRY
    if(!SysParams::checkChemParameters(ChemData))
        exit(EXIT_FAILURE);
#endif
#ifdef MECHANICS
    if(!SysParams::checkMechParameters(MTypes))
        exit(EXIT_FAILURE);
#endif
#ifdef DYNAMICRATES
    if(!SysParams::checkDyRateParameters(DRTypes))
        exit(EXIT_FAILURE);
#endif
    
    cout << "Done." << endl;
    
    //Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
    vector<vector<vector<double>>> filamentData;
    cout << "---" << endl;
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
    cout << "Done. " << filamentData.size()
         << " filaments created." << endl;
    
    //First update of system
    updatePositions();
    updateNeighborLists();
}

void Controller::updatePositions() {
    
    //update all moveables
    
    //Beads
    for(auto &f : *FilamentDB::instance()) {
        for (auto cylinder : f->getCylinderVector()){
            cylinder->getFirstBead()->updatePosition();
        }
        //update last bead
        f->getCylinderVector().back()->
        getSecondBead()->updatePosition();
    }
    
    //Cylinders
    for(auto &f : *FilamentDB::instance()) {
        for (auto cylinder : f->getCylinderVector()){
            cylinder->updatePosition();
        }
    }
    //Linkers
    for(auto &l : *LinkerDB::instance())
        l->updatePosition();
    
    //Motors
    for(auto &m : *MotorGhostDB::instance())
        m->updatePosition();
    
    //Branchers
    for(auto &b : *BranchingPointDB::instance())
        b->updatePosition();
}

#ifdef DYNAMICRATES
void Controller::updateReactionRates() {
    /// update all reactables
    
    //Boundary cylinders
    for(auto &c : _subSystem->getBoundaryCylinders())
        c->updateReactionRates();
    
    //Linkers
    for(auto &l : *LinkerDB::instance())
        l->updateReactionRates();
    
    //Motors
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
    
#ifdef CHEMISTRY
    double tauLastSnapshot = 0;
    double oldTau = 0;
#endif
    
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
    
    ///Print initial configuration
    for(auto o: _outputs) o->print(0);
    
    cout << "---" << endl;
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
            if(!_cController->run(_numChemSteps)) {
                for(auto o: _outputs) o->print(i + _numChemSteps);
                break;
            }
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

#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(i % _numStepsPerNeighbor == 0 && i != 0)
                updateNeighborLists(true);
            
            i += _numChemSteps;
            oldTau = tau();
        }
#endif
    }
    
    else {
#ifdef CHEMISTRY
        for(int i = 0; i < _numTotalSteps; i+=_numChemSteps) {
            //run ccontroller
            if(!_cController->run(_numChemSteps)) {
                for(auto o: _outputs) o->print(i + _numChemSteps);
                break;
            }
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

#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(i % _numStepsPerNeighbor == 0 && i != 0)
                updateNeighborLists(true);
            
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



