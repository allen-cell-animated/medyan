
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

#include "Parser.h"
#include "Output.h"
#include "SubSystem.h"
#include "BoundaryImpl.h"
#include "FilamentInitializer.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Controller::Controller(SubSystem* s) : _subSystem(s) {
    
    //init subsystem
    _subSystem = new SubSystem();
    
    //init controllers
    _mController   = new MController(_subSystem);
    _cController   = new CController(_subSystem);
    _gController   = new GController();
    _drController  = new DRController();
    
    //set Trackable's subsystem ptr
    Trackable::_subSystem = _subSystem;
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
    _inputDirectory  = inputDirectory;
    _outputDirectory = outputDirectory;
    
    //Parse input, get parameters
    SystemParser p(inputFile);
    
    //read output types
    OutputTypes oTypes = p.readOutputTypes();
    
    //snapshot type output
    cout << endl;
    
    string snapName   = _outputDirectory + "snapshot.traj";
    string birthName  = _outputDirectory + "birthtimes.traj";
    string forceName  = _outputDirectory + "forces.traj";
    string stressName = _outputDirectory + "stresses.traj";
    
    if(oTypes.basicSnapshot)
        _outputs.push_back(new BasicSnapshot(snapName));
    if(oTypes.birthTimes)
        _outputs.push_back(new BirthTimes(birthName));
    if(oTypes.forces)
        _outputs.push_back(new Forces(forceName));
    if(oTypes.stresses)
        _outputs.push_back(new Stresses(stressName));
    
    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters())
        exit(EXIT_FAILURE);
    
    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "---" << endl;
    cout << "Initializing geometry...";
    CompartmentGrid* grid = _gController->initializeGrid();
    _subSystem->setCompartmentGrid(grid);
    cout << "Done." << endl;
    
    //Always read boundary type
    auto BTypes = p.readBoundaryType();
    p.readBoundParams();
    
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
    
    //Initialize boundary
    cout << "---" << endl;
    cout << "Initializing boundary...";
    if(BTypes.boundaryShape == "CUBIC") {
        _subSystem->addBoundary(new BoundaryCubic(_subSystem));
    }
    else if(BTypes.boundaryShape == "SPHERICAL") {
        _subSystem->addBoundary(
        new BoundarySpherical(_subSystem, SysParams::Boundaries().diameter));
    }
    else if(BTypes.boundaryShape == "CAPSULE") {
        _subSystem->addBoundary(
        new BoundaryCapsule(_subSystem, SysParams::Boundaries().diameter));
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
    _numChemSteps = CAlgorithm.numChemSteps;
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
    cout << "Checking cross-parameter consistency...";
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
    for (auto it: filamentData) {
        
        double d = twoPointDistance(it[0], it[1]);
        vector<double> tau = twoPointDirection(it[0], it[1]);

        int numSegment = d / SysParams::Geometry().cylinderSize;

        // check how many segments can fit between end-to-end of the filament
        if (numSegment == 0)
            _subSystem->addTrackable<Filament>(_subSystem, it, 2);
        else
            _subSystem->addTrackable<Filament>(_subSystem, it, numSegment + 1);
    }
    cout << "Done. " << filamentData.size() << " filaments created." << endl;
}

void Controller::updatePositions() {
    
    //NEED TO UPDATE CYLINDERS FIRST
    for(auto c : Cylinder::getCylinders())
        c->updatePosition();
    
    //update all other moveables
    for(auto m : _subSystem->getMovables())
        m->updatePosition();
}

#ifdef DYNAMICRATES
void Controller::updateReactionRates() {
    /// update all reactables
    for(auto r : _subSystem->getReactables())
        r->updateReactionRates();
}
#endif

void Controller::updateNeighborLists() {
    
    //Full reset of neighbor lists
    _subSystem->resetNeighborLists();
    
#ifdef CHEMISTRY
    _subSystem->updateBindingManagers();
#endif
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
    
    //reupdate positions and neighbor lists
    updatePositions();
    updateNeighborLists();
    
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
            if(i % _numStepsPerNeighbor == 0)
                updateNeighborLists();
            
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
            if(i % _numStepsPerNeighbor == 0)
                updateNeighborLists();
            
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



