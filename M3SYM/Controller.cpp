
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
#include "Boundary.h"
#include "CompartmentGrid.h"

#include "FilamentInitializer.h"
#include "BubbleInitializer.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "MTOC.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Controller::Controller(SubSystem* s) : _subSystem(s) {
    
    //init subsystem
    _subSystem = new SubSystem();
    
    //init controllers
    _mController   = new MController(_subSystem);
    _cController   = new CController(_subSystem);
    _gController   = new GController(_subSystem);
    _drController  = new DRController();
    
    //set Trackable's subsystem ptr
    Trackable::_subSystem = _subSystem;
}

void Controller::initialize(string inputFile,
                            string inputDirectory,
                            string outputDirectory) {
    
    ///Seed simple rand generation for small tasks
    srand(time(NULL));
    
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
    string tensionName = _outputDirectory + "tensions.traj";
    
    if(oTypes.basicSnapshot) _outputs.push_back(new BasicSnapshot(snapName));
    if(oTypes.birthTimes)    _outputs.push_back(new BirthTimes(birthName));
    if(oTypes.forces)        _outputs.push_back(new Forces(forceName));
    if(oTypes.tensions)      _outputs.push_back(new Tensions(tensionName));
    
    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters()) exit(EXIT_FAILURE);
    
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
    //Initialize boundary
    cout << "---" << endl;
    cout << "Initializing boundary...";
    
    auto BTypes = p.readBoundaryType();
    p.readBoundParams();
    
    //initialize
    _gController->initializeBoundary(BTypes);
    cout << "Done." <<endl;
    
#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController->setActiveCompartments();
    
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
    
    //Set up chemistry output if any
    string chemsnapname = _outputDirectory + "chemistry.traj";
    if(oTypes.chemistry)
        _outputs.push_back(new Chemistry(chemsnapname, ChemData,
                                         _subSystem->getCompartmentGrid()));
#endif
    
#ifdef DYNAMICRATES
    cout << "---" << endl;
    cout << "Initializing dynamic rates...";
    //read dynamic rate parameters
    p.readDyRateParams();
    
    //read dynamic rate types
    DynamicRateType DRTypes = p.readDynamicRateType();
    
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
    
    //setup initial network configuration
    setupInitialNetwork(p);
    
    //setup special structures
    setupSpecialStructures(p);
}

void Controller::setupInitialNetwork(SystemParser& p) {
    
    //Read bubble setup, parse bubble input file if needed
    BubbleSetup BSetup = p.readBubbleSetup();
    BubbleData bubbles;
    
    cout << "---" << endl;
    cout << "Initializing bubbles...";
    
    if(BSetup.inputFile != "") {
        BubbleParser bp(_inputDirectory + BSetup.inputFile);
        bubbles = bp.readBubbles();
    }
    //add other bubbles if specified
    BubbleInitializer* bInit = new RandomBubbleDist();
    
    auto bubblesGen = bInit->createBubbles(_subSystem->getBoundary(),
                                           BSetup.numBubbles,
                                           BSetup.bubbleType);
    bubbles.insert(bubbles.end(), bubblesGen.begin(), bubblesGen.end());
    delete bInit;
    
    //add bubbles
    for (auto it: bubbles) {
        
        auto coord = get<1>(it);
        auto type = get<0>(it);
        
        if(type >= SysParams::Mechanics().numBubbleTypes) {
            cout << "Bubble data specified contains an "
                 <<"invalid bubble type. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        _subSystem->addTrackable<Bubble>(_subSystem, coord, type);
    }
    cout << "Done. " << bubbles.size() << " bubbles created." << endl;
    
    //Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
    FilamentData filaments;
    
    cout << "---" << endl;
    cout << "Initializing filaments...";
    
    if(FSetup.inputFile != "") {
        FilamentParser fp(_inputDirectory + FSetup.inputFile);
        filaments = fp.readFilaments();
    }
    
    //add other filaments if specified
    FilamentInitializer* fInit = new RandomFilamentDist();
    
    auto filamentsGen = fInit->createFilaments(_subSystem->getBoundary(),
                                               FSetup.numFilaments,
                                               FSetup.filamentType,
                                               FSetup.filamentLength);
    
    filaments.insert(filaments.end(), filamentsGen.begin(), filamentsGen.end());
    delete fInit;
    
    //add filaments
    for (auto it: filaments) {
        
        auto coord1 = get<1>(it);
        auto coord2 = get<2>(it);
        auto type = get<0>(it);
        
        if(type >= SysParams::Chemistry().numFilaments) {
            cout << "Filament data specified contains an "
                 <<"invalid filament type. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        vector<vector<double>> coords = {coord1, coord2};
        
        double d = twoPointDistance(coord1, coord2);
        vector<double> tau = twoPointDirection(coord1, coord2);
        
        int numSegment = d / SysParams::Geometry().cylinderSize[type];
        
        // check how many segments can fit between end-to-end of the filament
        if (numSegment == 0)
            _subSystem->addTrackable<Filament>(_subSystem, type, coords, 2);
        else
            _subSystem->addTrackable<Filament>(_subSystem, type, coords, numSegment + 1);
    }
    cout << "Done. " << filaments.size() << " filaments created." << endl;
}

void Controller::setupSpecialStructures(SystemParser& p) {
    
    cout << "---" << endl;
    cout << "Setting up special structures...";
    
    SpecialSetupType SType = p.readSpecialSetupType();
    
    //set up a MTOC if desired
    //For now, uses 20 filaments
    if(SType.mtoc) {
        
        MTOC* mtoc = _subSystem->addTrackable<MTOC>();
        
        //create the bubble in top part of grid, centered in x,y
        double bcoordx = GController::getSize()[0] / 2;
        double bcoordy = GController::getSize()[1] / 2;
        double bcoordz = GController::getSize()[2] * 5 / 6;
        
        vector<double> bcoords = {bcoordx, bcoordy, bcoordz};
        Bubble* b = _subSystem->addTrackable<Bubble>(_subSystem, bcoords, SType.mtocBubbleType);

        mtoc->setBubble(b);
        
        FilamentInitializer *init = new MTOCFilamentDist(bcoords,
        SysParams::Mechanics().BubbleRadius[SType.mtocBubbleType]);
        
        auto filaments = init->createFilaments(_subSystem->getBoundary(),
                                               SType.mtocNumFilaments,
                                               SType.mtocFilamentType,
                                               SType.mtocFilamentLength);
        //add filaments
        for (auto it: filaments) {
            
            auto coord1 = get<1>(it);
            auto coord2 = get<2>(it);
            
            vector<vector<double>> coords = {coord1, coord2};
            
            double d = twoPointDistance(coord1, coord2);
            vector<double> tau = twoPointDirection(coord1, coord2);
            
            int numSegment = d / SysParams::Geometry().cylinderSize[SType.mtocFilamentType];
            
            // check how many segments can fit between end-to-end of the filament
            Filament *f = _subSystem->addTrackable<Filament>(_subSystem, SType.mtocFilamentType,
                                                             coords, numSegment + 1, "ARC");
            
            mtoc->addFilament(f);
        }
    }
    cout << "Done." << endl;
}

void Controller::moveBoundary(double deltaTau) {
    
    //calculate distance to move
    double dist = SysParams::Boundaries().moveSpeed * deltaTau;
    
    //move it
    if(tau() >= SysParams::Boundaries().moveStartTime &&
       tau() <= SysParams::Boundaries().moveEndTime)
        _subSystem->getBoundary()->move(dist);
    
    //activate, deactivate necessary compartments
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
        
        if(_subSystem->getBoundary()->within(C)) {
            
            if(C->isActivated()) continue;
            else _cController->activate(C);
        }
        else {
            if(!C->isActivated()) continue;
            else _cController->deactivate(C);
        }
    }
}

void Controller::updatePositions() {
    
    //NEED TO UPDATE CYLINDERS FIRST
    for(auto c : Cylinder::getCylinders()) c->updatePosition();
    
    //update all other moveables
    for(auto m : _subSystem->getMovables()) m->updatePosition();
}

#ifdef DYNAMICRATES
void Controller::updateReactionRates() {
    /// update all reactables
    for(auto r : _subSystem->getReactables()) r->updateReactionRates();
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
    
    cout << "---" << endl;
    cout << "Performing an initial minimization..." << endl;
    
    //perform first minimization
#ifdef MECHANICS
    _mController->run(false);
    
    //reupdate positions and neighbor lists
    updatePositions();
    updateNeighborLists();
    
#ifdef DYNAMICRATES
    updateReactionRates();
#endif
    
#endif
    ///Print initial configuration
    for(auto o: _outputs) o->print(0);
    
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
            i += _numChemSteps;
            
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
            
            //move the boundary
            moveBoundary(tau() - oldTau);
            
            oldTau = tau();
            
            _numTotalSteps = i;
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
            
            //move the boundary
            moveBoundary(tau() - oldTau);
            
            oldTau = tau();
        }
#endif
    }
    //print last snapshots
    for(auto o: _outputs) o->print(_numTotalSteps);
    
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
}



