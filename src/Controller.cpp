
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
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
#include "Structure/BoundaryElementImpl.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "MTOC.h"
#include "AFM.h"
#include "ChemManager.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "MController.h"
#include "Cylinder.h"
#include <unordered_map>
#include     <tuple>
#include <vector>
#include <algorithm>
#include "ChemManager.h"

#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Util/Io/Log.hpp"
#include "Util/Profiler.hpp"

using namespace mathfunc;

Controller::Controller() :
    _mController(&_subSystem),
    _cController(&_subSystem),
    _gController(&_subSystem) {

    //set Trackable's subsystem ptr
    Trackable::_subSystem = &_subSystem;
}

void Controller::initialize(string inputFile,
                            string inputDirectory,
                            string outputDirectory,
                            int numThreads) {

    // Notice: numThreads is not used currently.

    SysParams::INITIALIZEDSTATUS = false;
    //general check of macros
#if defined(DYNAMICRATES) && (!defined(CHEMISTRY) || !defined(MECHANICS))
    LOG(FATAL) << "If dynamic rates is turned on, chemistry and mechanics must be "
         << "defined. Please set these compilation macros and try again. Exiting.";
    exit(EXIT_FAILURE);
#endif

    //init input directory
    _inputDirectory  = inputDirectory + "/";
    _outputDirectory = outputDirectory + "/";

    //Parse input, get parameters
    _inputFile = inputFile;
    SystemParser p(inputFile);

    //snapshot type output
    cout << endl;

    //trajectory-style data
    _outputs.push_back(make_unique<BasicSnapshot>(_outputDirectory + "snapshot.traj", &_subSystem));
    _outputs.push_back(make_unique<BirthTimes>(_outputDirectory + "birthtimes.traj", &_subSystem));
    _outputs.push_back(make_unique<Forces>(_outputDirectory + "forces.traj", &_subSystem));
    _outputs.push_back(make_unique<Tensions>(_outputDirectory + "tensions.traj", &_subSystem));

    _outputs.push_back(make_unique<PlusEnd>(_outputDirectory + "plusend.traj", &_subSystem));
    //ReactionOut should be the last one in the output list
    //Otherwise incorrect deltaMinusEnd or deltaPlusEnd values may be genetrated.
    _outputs.push_back(make_unique<ReactionOut>(_outputDirectory + "monomers.traj", &_subSystem));
    //add br force out and local diffussing species concentration
    _outputs.push_back(make_unique<BRForces>(_outputDirectory + "repulsion.traj", &_subSystem));
    //_outputs.push_back(make_unique<PinForces>(_outputDirectory + "pinforce.traj", &_subSystem));

    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters()) exit(EXIT_FAILURE);

    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "---" << endl;
    LOG(STEP) << "Initializing geometry...";
    _gController.initializeGrid();
    LOG(INFO) << "Done.";

    //Initialize boundary
    cout << "---" << endl;
    LOG(STEP) << "Initializing boundary...";

    auto BTypes = p.readBoundaryType();
    p.readBoundParams();

    //initialize
    _gController.initializeBoundary(BTypes);
    LOG(INFO) << "Done.";

#ifdef MECHANICS
    //read algorithm and types
    auto MTypes = p.readMechanicsFFType();
    auto MAlgorithm = p.readMechanicsAlgorithm();

    //read const parameters
    p.readMechParams();

    //Initialize Mechanical controller
    cout << "---" << endl;
    LOG(STEP) << "Initializing mechanics...";
    _mController.initialize(MTypes, MAlgorithm);
    LOG(INFO) << "Done.";

#endif

#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController.setActiveCompartments();

    if(_subSystem.getBoundary()->getShape() == BoundaryShape::Cylinder){
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()){
            C->computeSlicedVolumeArea(Compartment::SliceMethod::cylinderBoundary);
        }
    }
    else{
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()){
            C->computeNonSlicedVolumeArea();
        }
    }
    //Calculate surface area and volume for reaction rate scaling


    //read parameters
    p.readChemParams();

    //Initialize chemical controller
    cout << "---" << endl;
    LOG(STEP) << "Initializing chemistry...";
    //read algorithm
    auto CAlgorithm = p.readChemistryAlgorithm();
    auto CSetup = p.readChemistrySetup();
    _cAlgorithm=CAlgorithm;
    //run time for sim
    _runTime = CAlgorithm.runTime;

    //freq of snapshots, minimizations, neighborlist updates
    _snapshotTime = CAlgorithm.snapshotTime;
    _minimizationTime = CAlgorithm.minimizationTime;
    _neighborListTime = CAlgorithm.neighborListTime;
    _datadumpTime = CAlgorithm.datadumpTime;

    //if run time was not set, look for runsteps parameters
    _runSteps = CAlgorithm.runSteps;
    _snapshotSteps = CAlgorithm.snapshotSteps;
    _minimizationSteps = CAlgorithm.minimizationSteps;
    _neighborListSteps = CAlgorithm.neighborListSteps;

    ChemistryData ChemData;

    if(CSetup.inputFile != "") {
        ChemistryParser cp(_inputDirectory + CSetup.inputFile);
        ChemData = cp.readChemistryInput();
        _chemData=ChemData;
    }
    else {
        LOG(FATAL) << "Need to specify a chemical input file. Exiting.";
        exit(EXIT_FAILURE);
    }

#ifdef CHEMISTRY
    SysParams::addChemParameters(ChemData);

    if(!SysParams::checkChemParameters(ChemData))
        exit(EXIT_FAILURE);
#endif

    // create the dissiption tracking object
    _dt = new DissipationTracker(&_mController);
    _cController.initialize(CAlgorithm.algorithm, ChemData, _dt);
    LOG(INFO) << "Done.";

    //Set up chemistry output if any
    string chemsnapname = _outputDirectory + "chemistry.traj";
    _outputs.push_back(make_unique<Chemistry>(chemsnapname, &_subSystem, ChemData,
                                     _subSystem.getCompartmentGrid()));

    ChemSim* _cs = _cController.getCS();
	ForceFieldManager* _ffm = _mController.getForceFieldManager();

    string concenname = _outputDirectory + "concentration.traj";

    _outputs.push_back(make_unique<Concentrations>(concenname, &_subSystem, ChemData));

    if(SysParams::CParams.dissTracking){
        //Set up dissipation output if dissipation tracking is enabled
        string disssnapname = _outputDirectory + "dissipation.traj";
        _outputs.push_back(make_unique<Dissipation>(disssnapname, &_subSystem, _cs));

        //Set up HRCD output if dissipation tracking is enabled
        string hrcdsnapname = _outputDirectory + "HRCD.traj";

        _outputs.push_back(make_unique<HRCD>(hrcdsnapname, &_subSystem, _cs));

        //Set up HRMD output if dissipation tracking is enabled
        string hrmdsnapname = _outputDirectory + "HRMD.traj";
        _outputs.push_back(make_unique<HRMD>(hrmdsnapname, &_subSystem, _cs));

    }

    if(SysParams::CParams.eventTracking){
        //Set up MotorWalkingEvents if event tracking is enabled
        string motorwalkingevents = _outputDirectory + "motorwalkingevents.traj";
        _outputs.push_back(make_unique<MotorWalkingEvents>(motorwalkingevents, &_subSystem, _cs));

        //Set up motorunbindingevents if event tracking is enabled
        string motorunbindingevents = _outputDirectory + "motorunbindingevents.traj";
        _outputs.push_back(make_unique<MotorUnbindingEvents>(motorunbindingevents, &_subSystem, _cs));

        //Set up LinkerUnbindingEvents if event tracking is enabled
        string linkerunbindingevents = _outputDirectory + "linkerunbindingevents.traj";
        _outputs.push_back(make_unique<LinkerUnbindingEvents>(linkerunbindingevents, &_subSystem, _cs));

        //Set up LinkerBindingEvents if event tracking is enabled
        string linkerbindingevents = _outputDirectory + "linkerbindingevents.traj";
        _outputs.push_back(make_unique<LinkerBindingEvents>(linkerbindingevents, &_subSystem, _cs));
    }

    if(SysParams::MParams.hessTracking){
        //Set up HessianMatrix if hessiantracking is enabled
        string hessianmatrix = _outputDirectory + "hessianmatrix.traj";
        _outputs.push_back(make_unique<HessianMatrix>(hessianmatrix, &_subSystem, _ffm));

        //Set up HessianSpectra if hessiantracking is enabled
        string hessianspectra = _outputDirectory + "hessianspectra.traj";
        _outputs.push_back(make_unique<HessianSpectra>(hessianspectra, &_subSystem, _ffm));

    }

    //Set up CMGraph output
    string cmgraphsnapname = _outputDirectory + "CMGraph.traj";
    _outputs.push_back(make_unique<CMGraph>(cmgraphsnapname, &_subSystem));
    
    //Set up TMGraph output
    string tmgraphsnapname = _outputDirectory + "TMGraph.traj";
    _outputs.push_back(make_unique<TMGraph>(tmgraphsnapname, &_subSystem));


    //Set up datadump output if any
    string datadumpname = _outputDirectory + "datadump.traj";
    _outputdump.push_back(make_unique<Datadump>(datadumpname, &_subSystem, ChemData));

//    string twofilamentname = _outputDirectory + "twofilament.traj";
//    _outputs.push_back(make_unique<TwoFilament>(twofilamentname, &_subSystem, ChemData));

//    //Set up Turnover output if any
//    string turnover = _outputDirectory + "Turnover.traj";
//    _outputs.push_back(make_unique<FilamentTurnoverTimes>(turnover, &_subSystem));


#endif

#ifdef DYNAMICRATES
    cout << "---" << endl;
    LOG(STEP) << "Initializing dynamic rates...";
    //read dynamic rate parameters
    p.readDyRateParams();

    //read dynamic rate types
    DynamicRateType DRTypes = p.readDynamicRateType();

    //init controller
    _drController.initialize(DRTypes);
    LOG(INFO) << "Done.";

#endif

    //Check consistency of all chemistry and mechanics parameters
    cout << "---" << endl;
    LOG(STEP) << "Checking cross-parameter consistency...";
    //Chemistry is checked in advance
#ifdef MECHANICS
    if(!SysParams::checkMechParameters(MTypes))
        exit(EXIT_FAILURE);
#endif
#ifdef DYNAMICRATES
    if(!SysParams::checkDyRateParameters(DRTypes))
        exit(EXIT_FAILURE);
#endif

    LOG(INFO) << "Done.";

    //setup initial network configuration
    setupInitialNetwork(p);

    //setup special structures
    p.readSpecialParams();
    setupSpecialStructures(p);

    SysParams::INITIALIZEDSTATUS = true;
}

void Controller::setupInitialNetwork(SystemParser& p) {

    //Read bubble setup, parse bubble input file if needed
    BubbleSetup BSetup = p.readBubbleSetup();
    BubbleData bubbles;

    cout << "---" << endl;
    cout << "Initializing bubbles...";

    if (BSetup.inputFile != "") {
        BubbleParser bp(_inputDirectory + BSetup.inputFile);
        bubbles = bp.readBubbles();
    }
    //add other bubbles if specified
    BubbleInitializer *bInit = new RandomBubbleDist();

    auto bubblesGen = bInit->createBubbles(_subSystem.getBoundary(),
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
        _subSystem.addTrackable<Bubble>(&_subSystem, coord, type);
    }
    cout << "Done. " << bubbles.size() << " bubbles created." << endl;

    //Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
//    FilamentData filaments;

    cout << "---" << endl;
//    HybridBindingSearchManager::setdOut();
    cout << "Initializing filaments...";

    if (SysParams::RUNSTATE == true) {
        if (FSetup.inputFile != "") {
            FilamentParser fp(_inputDirectory + FSetup.inputFile);
            filaments = fp.readFilaments();
        }
        fil = get<0>(filaments);
        //add other filaments if specified
        FilamentInitializer *fInit = new RandomFilamentDist();

        auto filamentsGen = fInit->createFilaments(_subSystem.getBoundary(),
                                                    FSetup.numFilaments,
                                                    FSetup.filamentType,
                                                    FSetup.filamentLength);
        auto filGen = get<0>(filamentsGen);
        fil.insert(fil.end(), filGen.begin(), filGen.end());
        delete fInit;

        //add filaments

        for (auto it: fil) {

            auto coord1 = get<1>(it);
            auto coord2 = get<2>(it);
            auto type = get<0>(it);

            if (type >= SysParams::Chemistry().numFilaments) {
                cout << "Filament data specified contains an "
                        << "invalid filament type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            vector<vector<floatingpoint>> coords = {coord1, coord2};

            if (coord2.size() == 3) {

                floatingpoint d = twoPointDistance(coord1, coord2);
                vector<floatingpoint> tau = twoPointDirection(coord1, coord2);
                int numSegment = static_cast<int>(std::round(
                        d / SysParams::Geometry().cylinderSize[type]));

                // check how many segments can fit between end-to-end of the filament
                if (numSegment == 0)
                    _subSystem.addTrackable<Filament>(&_subSystem, type, coords, 2,
                                                        FSetup.projectionType);
                else
                    _subSystem.addTrackable<Filament>(&_subSystem, type, coords,
                                                        numSegment + 1,
                                                        FSetup.projectionType);
            } else if (coord2.size() > 3) {
                int numSegment = coord2.size() / 3;
                vector<vector<floatingpoint>> coords;
                coords.push_back(coord1);
                for (int id = 0; id < numSegment; id++)
                    coords.push_back(
                            {coord2[id * 3], coord2[id * 3 + 1], coord2[id * 3 + 2]});

                if (numSegment == 0)
                    _subSystem.addTrackable<Filament>(&_subSystem, type, coords, 2,
                                                        FSetup.projectionType);
                else
                    _subSystem.addTrackable<Filament>(&_subSystem, type, coords,
                                                        numSegment + 1,
                                                        FSetup.projectionType);
            }
        }
        cout << "Done. " << fil.size() << " filaments created." << endl;
        cout << "Total cylinders " << Cylinder::getCylinders().size() << endl;
    }
    else{
        cout<<endl;
	    cout<<"RESTART PHASE BEINGS."<<endl;
        //Create the restart pointer
        const string inputfileName = _inputDirectory + FSetup.inputFile;
        _restart = new Restart(&_subSystem, _chemData, inputfileName);
        //read set up.
        _restart->readNetworkSetup();
        _restart->setupInitialNetwork();
    }
}

void Controller::setupSpecialStructures(SystemParser& p) {

    cout << "---" << endl;
    cout << "Setting up special structures...";

    SpecialSetupType SType = p.readSpecialSetupType();

    //set up a MTOC if desired

    //For now, uses 20 filaments
    if(SType.mtoc) {

        MTOC* mtoc = _subSystem.addTrackable<MTOC>();
        
        //set MTOC coordinates based on input
        floatingpoint bcoordx = SType.mtocInputCoordXYZ[0];
        floatingpoint bcoordy = SType.mtocInputCoordXYZ[1];
        floatingpoint bcoordz = SType.mtocInputCoordXYZ[2];


        vector<floatingpoint> bcoords = {bcoordx, bcoordy, bcoordz};
        Bubble* b = _subSystem.addTrackable<Bubble>(&_subSystem, bcoords, SType.mtocBubbleType);


        mtoc->setBubble(b);
        
        FilamentInitializer *init = new MTOCFilamentDist(bcoords,
                                                         SysParams::Mechanics().BubbleRadius[SType.mtocBubbleType]);


        auto filaments = init->createFilaments(_subSystem.getBoundary(),

                                               SType.mtocNumFilaments,
                                               SType.mtocFilamentType,
                                               SType.mtocFilamentLength);
        //add filaments
        filamentData fil=get<0>(filaments);
        for (auto it: fil) {
            
            auto coord1 = get<1>(it);
            auto coord2 = get<2>(it);
            
            vector<vector<floatingpoint>> coords = {coord1, coord2};
            
            floatingpoint d = twoPointDistance(coord1, coord2);
            vector<floatingpoint> tau = twoPointDirection(coord1, coord2);
            
            int numSegment = d / SysParams::Geometry().cylinderSize[SType.mtocFilamentType];
            
            // check how many segments can fit between end-to-end of the filament
            Filament *f = _subSystem.addTrackable<Filament>(&_subSystem, SType.mtocFilamentType,
                                                             coords, numSegment + 1, "ARC");
            
            mtoc->addFilament(f);
            
        }
        cout << "MTOC is set." << endl;
        
    }
    else if(SType.afm) {

        AFM* afm = _subSystem.addTrackable<AFM>();

        //create a bubble in top part of grid, centered in x,y
        floatingpoint bcoordx = GController::getSize()[0] / 2;
        floatingpoint bcoordy = GController::getSize()[1] / 2;
        //set up the height of the AFM bubble
        floatingpoint bcoordz = 1250;

        vector<floatingpoint> bcoords = {bcoordx, bcoordy, bcoordz};
        Bubble* b = _subSystem.addTrackable<Bubble>(&_subSystem, bcoords, SType.afmBubbleType);

        PlaneBoundaryElement* afmpbe = _subSystem.addTrackable<PlaneBoundaryElement>(bcoords, vector<floatingpoint>{0,0,-1}, SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength);

        afm->setBubble(b);
        afm->setPlaneBoundaryElement(afmpbe);

        FilamentInitializer *init = new AFMFilamentDist(bcoords, SysParams::Mechanics().BubbleRadius[SType.afmBubbleType]);

        auto filaments = init->createFilaments(_subSystem.getBoundary(),
                                               SType.afmNumFilaments,
                                               SType.afmFilamentType,
                                               SType.afmFilamentLength);
        //add filaments
        filamentData fil=get<0>(filaments);
        for (auto it: fil) {

            auto coord1 = get<1>(it);
            auto coord2 = get<2>(it);

            vector<vector<floatingpoint>> coords = {coord1, coord2};

            floatingpoint d = twoPointDistance(coord1, coord2);
            vector<floatingpoint> tau = twoPointDirection(coord1, coord2);

            int numSegment = static_cast<int>(std::round(d / SysParams::Geometry().cylinderSize[SType.afmFilamentType]));

            // check how many segments can fit between end-to-end of the filament



            Filament *f = _subSystem.addTrackable<Filament>(&_subSystem, SType.afmFilamentType, coords, numSegment + 1, "ARC");

            afm->addFilament(f);
        }
        cout << "AFM is set." << endl;
    }
    cout << "Done." << endl;
}

void Controller::activatedeactivateComp(){

    if(SysParams::Boundaries().transfershareaxis>=0){
        fCompmap.clear();
        bCompmap.clear();
        activatecompartments.clear();
        ControlfrontbackEndComp();
//            std::cout<<fCompmap.size()<<" "<<bCompmap.size()<<" "<<activatecompartments.size()<<endl;
        for(auto it=activatecompartments.begin();it!=activatecompartments.end();it++)
        {
            if(!(*it)->isActivated())
                _cController.activate(*it);
        }
        //deactivate compartments starting from the right extreme
        for (std::multimap<int,Compartment*>::reverse_iterator it=fCompmap.rbegin(); it!=fCompmap.rend(); ++it)
            _cController.deactivate(it->second);
        //deactivate compartments starting from the left extreme
        for (std::multimap<int,Compartment*>::iterator it=bCompmap.begin(); it!=bCompmap.end(); ++it)
            _cController.deactivate(it->second);
        fCompmap.clear();
        bCompmap.clear();

//            std::cout<<"Printing diffusing actin copy numbers."<<endl;
//
//            for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
//                for(auto sd : _chemData.speciesDiffusing) {
//                    string name = get<0>(sd);
//                    if(name.find("AD") != string::npos){
//                        auto s = C->findSpeciesByName(name);
//                        auto copyNum = s->getN();
//
//                        std::cout <<C->coordinates()[0]<<" "<<copyNum<<" ";
//                    }
//                }
//            }
//            std::cout<<endl;
    }
}

void Controller::ControlfrontbackEndComp(){
    Compartment* maxcomp=NULL;
    Compartment* mincomp=NULL;
    int tsaxis = SysParams::Boundaries().transfershareaxis;
    int planestomove = SysParams::Boundaries().planestomove;
    bool maxcompstate = false;
    bool mincompstate = false;
    if(planestomove == 2 || planestomove == 0) maxcompstate = true;
    if(planestomove == 2 || planestomove == 1) mincompstate = true;
    floatingpoint systemspan = 0.0;
    floatingpoint cmpsize = 0.0;
    if(tsaxis == 0)
    {systemspan = SysParams::Geometry().NX * SysParams::Geometry()
                                                              .compartmentSizeX;
    cmpsize = SysParams::Geometry().compartmentSizeX;}
    else if(tsaxis == 1) {
        systemspan = SysParams::Geometry().NY * SysParams::Geometry()
                .compartmentSizeY;
        cmpsize = SysParams::Geometry().compartmentSizeY;
    }
    else if(tsaxis == 2) {
        systemspan = SysParams::Geometry().NZ * SysParams::Geometry().compartmentSizeZ;
        cmpsize = SysParams::Geometry().compartmentSizeZ;
    }
    //copy vector to prevcopy
    bounds_prev[0] = bounds[0];bounds_prev[1] = bounds[1];
    bounds[0] = 0.0; bounds[1] =  systemspan;
    for(auto C : _subSystem.getCompartmentGrid()->getCompartments()){
        auto cyls=C->getCylinders();
        if(cyls.size()>0){
            //maxcomp refers to the compartment on the right extreme of reaction volume
            // that represents the current chemical boundary of the system.
            if(maxcompstate) {
                if (maxcomp == NULL)
                    maxcomp = C;
                else {
                    //get current maxcomp coordinates
                    auto mcoord = maxcomp->coordinates();
                    //get compartment coorinates
                    auto ccord = C->coordinates();
                    //compare to see if the compartment is further to the right of maxcomp.
                    if (mcoord[SysParams::Boundaries().transfershareaxis] <
                        ccord[SysParams::Boundaries()
                                .transfershareaxis])
                        maxcomp = C;
                }
            }
            //mincomp refers to the compartment on the left extreme of reaction volume
            // that represents the current chemical boundary of the system.
            if(mincompstate) {
                if (mincomp == NULL)
                    mincomp = C;
                else {
                    auto mcoord = mincomp->coordinates();
                    auto ccord = C->coordinates();
                    //compare to see if the compartment is further to the left of mincomp.
                    if (mcoord[SysParams::Boundaries().transfershareaxis] >
                        ccord[SysParams::Boundaries().transfershareaxis])
                        mincomp = C;
                }
            }
        }
    }

    if(maxcompstate) {
        std::cout<<"1 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        // front end is defined two compartments away from the current maxcomp.
        auto cmaxcomp = maxcomp->coordinates();
        //get the neighbor who is to the right of maxcomp.
        for (auto C:maxcomp->getNeighbours()) {
            auto cC = C->coordinates();
            if (cmaxcomp[SysParams::Boundaries().transfershareaxis] <
                cC[SysParams::Boundaries().transfershareaxis])
                maxcomp = C;
        }
        std::cout<<"2 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        cmaxcomp = maxcomp->coordinates();
        //get the neighbor who is to the right of maxcomp.
        for (auto C:maxcomp->getNeighbours()) {
            auto cC = C->coordinates();
            if (cmaxcomp[SysParams::Boundaries().transfershareaxis] <
                cC[SysParams::Boundaries().transfershareaxis])
                maxcomp = C;
        }
        std::cout<<"3 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        cmaxcomp = maxcomp->coordinates();
        assert((maxcomp != NULL) && "Non existent maxcomp. Exiting.");
        //Loop through compartments
        for (auto C : _subSystem.getCompartmentGrid()->getCompartments()) {
            auto cC = C->coordinates();
            //if compartment is to the right of maxcomp and activated, add to a vector to
            // deactivate later.
            if (cC[SysParams::Boundaries().transfershareaxis] >
                cmaxcomp[SysParams::Boundaries().transfershareaxis]) {
                if (C->isActivated())
                    fCompmap.insert(pair<int, Compartment *>(
                            cC[SysParams::Boundaries().transfershareaxis], C));
            }
                //if compartment is to the left of maxcomp and not active, add to a vector to
                // activate later.
            else {
                if (!(C->isActivated()))
                    activatecompartments.push_back(C);
            }
        }
        bounds[1] = maxcomp->coordinates()[SysParams::Boundaries().transfershareaxis] +
                cmpsize/2;
    }
    //back end is defined as the compartment that is two compartments to the left of
    // mincomp.
    if(mincompstate) {
        auto cmincomp = mincomp->coordinates();
        //get the neighbor who is to the left of mincomp.
        for (auto C:mincomp->getNeighbours()) {
            auto cC = C->coordinates();
            if (cmincomp[SysParams::Boundaries().transfershareaxis] >
                cC[SysParams::Boundaries().transfershareaxis])
                mincomp = C;
        }
        cmincomp = mincomp->coordinates();
        //get the neighbor who is to the left of mincomp.
        for (auto C:mincomp->getNeighbours()) {
            auto cC = C->coordinates();
            if (cmincomp[SysParams::Boundaries().transfershareaxis] >
                cC[SysParams::Boundaries().transfershareaxis])
                mincomp = C;
        }
        cmincomp = mincomp->coordinates();
        assert(mincomp != NULL && "Non existent mincomp. Exiting.");
        //Loop through compartments
        for (auto C : _subSystem.getCompartmentGrid()->getCompartments()) {
            auto cC = C->coordinates();
            //if compartment C is to the left of mincomp and was added to
            // activatecompartments vector, remove. If it is already active, add to a vector
            // to deactivate later.
            if (cC[SysParams::Boundaries().transfershareaxis] <
                cmincomp[SysParams::Boundaries().transfershareaxis]) {
                auto it = std::find(activatecompartments.begin(),
                                    activatecompartments.end(), C);
                if (it != activatecompartments.end())
                    activatecompartments.erase(it);
                if (C->isActivated()) {
                    bCompmap.insert(pair<int, Compartment *>(
                            cC[SysParams::Boundaries().transfershareaxis], C));
                }
            }
        }
        bounds[0] = mincomp->coordinates()[SysParams::Boundaries().transfershareaxis] -
                cmpsize/2;
    }
    //print the maximum (right boundary) and minimum (left boundary) compartment spans.
    std::cout<<"Maxbound "<<bounds[1]<<" Minbound "<<bounds[0]<<endl;
}

void Controller::moveBoundary(floatingpoint deltaTau) {
    //calculate distance to move
    floatingpoint dist = SysParams::Boundaries().moveSpeed * deltaTau;

    if(SysParams::Boundaries().transfershareaxis>=0){
        vector<floatingpoint> distvec= {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        if(SysParams::Boundaries().transfershareaxis == 0){
            distvec[0] = bounds[0] - bounds_prev[0];
            distvec[1] = bounds[1] - bounds_prev[1];
        }
        else if(SysParams::Boundaries().transfershareaxis == 1){
            distvec[2] = bounds[0] - bounds_prev[0];
            distvec[3] = bounds[1] - bounds_prev[1];
        }
        else if(SysParams::Boundaries().transfershareaxis == 2){
            distvec[4] = bounds[0] - bounds_prev[0];
            distvec[5] = bounds[1] - bounds_prev[1];
        }
        _subSystem.getBoundary()->move(distvec);
    }
        //deprecated not good to use.
    else if(abs(dist)>0){
        vector<floatingpoint> distvec = {dist, -dist, dist, -dist, dist, -dist};
        //move it
        if(tau() >= SysParams::Boundaries().moveStartTime &&
           tau() <= SysParams::Boundaries().moveEndTime)
            _subSystem.getBoundary()->move(distvec);

        //activate, deactivate necessary compartments
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()) {

            if(_subSystem.getBoundary()->within(C)) {

                if(C->isActivated()) continue;
                else _cController.activate(C);
            }
            else {
                if(!C->isActivated()) continue;
                else _cController.deactivate(C);
            }
        }
    }
    //calculate system volume.
    _subSystem.getBoundary()->volume();
}

void Controller::executeSpecialProtocols() {

    //making filaments static
    if(SysParams::Chemistry().makeFilamentsStatic &&
       SysParams::Chemistry().makeFilamentsStaticTime <= tau()) {

        //loop through all cylinders, passivate (de)polymerization
        for(auto c : Cylinder::getCylinders())
            c->getCCylinder()->passivatefilreactions();
    }

    //making linkers static
    if(SysParams::Chemistry().makeLinkersStatic &&
       SysParams::Chemistry().makeLinkersStaticTime <= tau()) {

        // loop through all linkers, passivate unbinding
        for(auto l: Linker::getLinkers())
            l->getCLinker()->getOffReaction()->passivateReaction();
    }


    if(SysParams::Mechanics().pinBoundaryFilaments &&
       tau() >= SysParams::Mechanics().pinTime) {

        pinBoundaryFilaments();
    }

    //Qin
    if(SysParams::Mechanics().pinLowerBoundaryFilaments &&
       tau() >= SysParams::Mechanics().pinTime) {

        pinLowerBoundaryFilaments();
    }
    
}

void Controller::updatePositions() {

	chrono::high_resolution_clock::time_point minsp, minep;
    //NEED TO UPDATE CYLINDERS FIRST
	minsp = chrono::high_resolution_clock::now();
    //Reset Cylinder update position state
    Cylinder::setpositionupdatedstate = false;
    for(auto c : Cylinder::getCylinders()) {
	    c->updatePosition();
    }
#ifdef OPTIMOUT
    cout<<"Cylinder position updated"<<endl;
#endif
    //Reset state to updated state
	Cylinder::setpositionupdatedstate = true;
	minep = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> compartment_update(minep - minsp);
    updatepositioncylinder += compartment_update.count();

    minsp = chrono::high_resolution_clock::now();
    //update all other moveables
//    for(auto m : _subSystem->getMovables()) m->updatePosition();
	int count = 0;
	for(auto m : Movable::getMovableList()) m->updatePosition();
    
    //update bubble
    if(SysParams::Chemistry().makeAFM) updateBubblePositions();


    minep = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> compartment_update2(minep - minsp);
    updatepositionmovable += compartment_update2.count();
}

void Controller::updateBubblePositions() {
    
    //update AFM bubble again based on time
    for(auto b : Bubble::getBubbles()) {
        if(b->isAFM()) b->updatePositionManually();
    }
    
    if(SysParams::Chemistry().makeRateDepend && tau() - tp > 1) {
        tp+=1;
        
        for(auto &filament : Filament::getFilaments()) {
            double deltaL;
            double numCyl = 0;
            for (auto cylinder : filament->getCylinderVector()){
                
                deltaL += cylinder->getMCylinder()->getLength() -
                cylinder->getMCylinder()->getEqLength();
                numCyl += 1;
            }
            
            //print last
            Cylinder* cylinder = filament->getCylinderVector().back();
            deltaL += cylinder->getMCylinder()->getLength() -
            cylinder->getMCylinder()->getEqLength();
            numCyl += 1;
            
            double k = cylinder->getMCylinder()->getStretchingConst();
            
            //if the filament tension is higher than threshold, regardless of sign
            if(k*deltaL/numCyl > SysParams::Chemistry().makeRateDependForce ||
               -k*deltaL/numCyl > SysParams::Chemistry().makeRateDependForce ){
                
                Cylinder* pCyl = filament->getCylinderVector().back();
                for(auto &r : pCyl->getCCylinder()->getInternalReactions()) {
                    if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
                        float newrate = 5 * SysParams::Chemistry().originalPolyPlusRate;
                        r->setBareRate(newrate);
                        r->recalcRateVolumeFactor();
                        r->updatePropensity();
                    }
                }
            }
            //else, set it back to orginal rate
            else{
                Cylinder* pCyl = filament->getCylinderVector().back();
                for(auto &r : pCyl->getCCylinder()->getInternalReactions()) {
                    if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
                        float newrate = SysParams::Chemistry().originalPolyPlusRate;
                        r->setBareRate(newrate);
                        r->recalcRateVolumeFactor();
                        r->updatePropensity();
                    }
                }
            }
        }
    }
}


#ifdef DYNAMICRATES
void Controller::updateReactionRates() {
    /// update all reactables
//    for(auto r : _subSystem->getReactables()) r->updateReactionRates();
	for(auto r : Reactable::getReactableList()){
		r->updateReactionRates();
		}
}
#endif

void Controller::updateNeighborLists() {
    #ifdef CROSSCHECK_CYLINDER
    if(HybridNeighborList::_crosscheckdumpFileNL.is_open())
        HybridNeighborList::_crosscheckdumpFileNL.close();
    string crosscheckNLname = _outputDirectory + "crosscheckNL.traj";
    HybridNeighborList::_crosscheckdumpFileNL.open(crosscheckNLname);
    #endif
    chrono::high_resolution_clock::time_point mins, mine;

    mins = chrono::high_resolution_clock::now();
    //Full reset of neighbor lists
    _subSystem.resetNeighborLists();
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runnl2(mine - mins);
    nl2time += elapsed_runnl2.count();
#ifdef CHEMISTRY
    mins = chrono::high_resolution_clock::now();
    _subSystem.updateBindingManagers();
#ifdef OPTIMOUT
	cout<<"updated BindingManagers"<<endl;
#endif
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runb(mine - mins);
    bmgrtime += elapsed_runb.count();
#endif
}

void Controller::resetCounters() {
    for(Filament* f : Filament::getFilaments()) f->resetCounters();
}

void Controller::pinBoundaryFilaments() {

    //if we've already added pinned filaments, return
    if(Bead::getPinnedBeads().size() != 0)
        return;

    //loop through beads, check if within pindistance
    for(auto b : Bead::getBeads()) {

        //pin only beads who are at the front of a plus end cylinder or back of a minus end cylinder
        Filament* f = (Filament*) b->getParent();
        Cylinder* plusEndC = f->getPlusEndCylinder();
        Cylinder* minusEndC = f->getMinusEndCylinder();

        if((plusEndC->getSecondBead() == b) ||
           (minusEndC->getFirstBead() == b)) {

            cout << _subSystem.getBoundary()->distance(b->vcoordinate()) << endl;
            cout << SysParams::Mechanics().pinDistance << endl;


            //if within dist to boundary, add
            if(_subSystem.getBoundary()->distance(b->vcoordinate()) < SysParams::Mechanics().pinDistance) {

                b->pinnedPosition = b->vcoordinate();
                b->addAsPinned();
            }
        }
    }
}

void Controller::pinLowerBoundaryFilaments() {

    //renew pinned filament list everytime

    //loop through beads, check if within pindistance
    for(auto b : Bead::getBeads()) {

        //pin all beads besides plus end and minus end cylinder
        Filament* f = (Filament*) b->getParent();
        Cylinder* plusEndC = f->getPlusEndCylinder();
        Cylinder* minusEndC = f->getMinusEndCylinder();

        if((plusEndC->getSecondBead() != b) ||
           (minusEndC->getFirstBead() != b)) {

            //cout << _subSystem->getBoundary()->lowerdistance(b->vcoordinate()) << endl;
            //cout << SysParams::Mechanics().pinDistance << endl;

            auto index = Rand::randfloatingpoint(0,1);
            //cout << index <<endl;
            //if within dist to boundary and index > 0.5, add
            if(_subSystem.getBoundary()->lowerdistance(b->vcoordinate()) < SysParams::Mechanics().pinDistance
               && index < SysParams::Mechanics().pinFraction && b->isPinned() == false) {
                //cout << index << endl;
                b->pinnedPosition = b->vcoordinate();
                b->addAsPinned();
            }
        }
    }
}

void Controller::run() {

#ifdef CHEMISTRY
    floatingpoint tauLastSnapshot = 0;
    floatingpoint tauLastMinimization = 0;
    floatingpoint tauLastNeighborList = 0;
    floatingpoint oldTau = 0;
    floatingpoint tauDatadump = 0;

    long stepsLastSnapshot = 0;
    long stepsLastMinimization = 0;
    long stepsLastNeighborList = 0;

    long totalSteps = 0;
#endif
    chrono::high_resolution_clock::time_point chk1, chk2, mins, mine;
    chk1 = chrono::high_resolution_clock::now();
//RESTART PHASE BEGINS
    if(SysParams::RUNSTATE==false){
//Step 2A. Turn off diffusion, passivate filament reactions and add reactions to heap.
        _restart->settorestartphase();
	    cout<<"Turned off Diffusion, and filament reactions."<<endl;
        cout<<"Bound species added to reaction heap."<<endl;

//Step 3. ############ RUN LINKER/MOTOR REACTIONS TO BIND BRANCHERS, LINKERS, MOTORS AT RESPECTIVE POSITIONS.#######
        cout<<"Number of reactions to be fired "<<_restart->getnumchemsteps()<<endl;
        _cController.runSteps(_restart->getnumchemsteps());
        cout<<"Reactions fired! Displaying number of reactions that are NOT fired in each"
              " compartment"
              ""<<endl;
//Step 4. Display the number of reactions yet to be fired. Should be zero.
        bool exitstatus = 0;
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto &Mgr:C->getFilamentBindingManagers()){
                int numsites = 0;
#ifdef NLORIGINAL
                numsites = Mgr->numBindingSites();
#else
                numsites = Mgr->numBindingSitesstencil();
#endif
                if(numsites == 0)
                    cout<< numsites<<" ";
                else{
                    cout<<endl;
                    LOG(ERROR)<<"Compartment ID "<<C->getId()<<" COORDS "
                                <<C->coordinates()[0] << " "
                                <<C->coordinates()[1] << " "
                                <<C->coordinates()[2] << endl;
                    LOG(ERROR)<<"Num binding sites "<<numsites<<endl;
                    string mgrname ="";
                    if(dynamic_cast<BranchingManager*>(Mgr.get()))
                        mgrname = " BRANCHING ";
                    else if (dynamic_cast<LinkerBindingManager*>(Mgr.get()))
                        mgrname = " LINKER ";
                    else
                        mgrname = " MOTOR ";
                    LOG(ERROR)<<"Printing "<<mgrname<<" binding sites that were not "
                                                      "chosen"<<endl;
                    #ifdef NLORIGINAL
                    Mgr->printbindingsites();
					#else
                	Mgr->printbindingsitesstencil();
					#endif
                	exitstatus = true;
                }
            }}
        cout<<endl;
        if(exitstatus) {
            cout << "Few reactions were not fired! Cannot restart this trajectory. "
                    "Exiting after printing diffusing species in each compartment..." <<
                    endl;

            cout<< "COMPARTMENT DATA: CMPID DIFFUSINGSPECIES COPYNUM"<<endl;
            for(auto cmp:_subSystem.getCompartmentGrid()->getCompartments()){
                cout <<cmp->getId()<<" ";
                for(auto sd : _chemData.speciesDiffusing) {
                    string name = get<0>(sd);
                    auto s = cmp->findSpeciesByName(name);
                    auto copyNum = s->getN();
                    cout <<name<<" "<<copyNum<<" ";
                }
                cout <<endl;
            }
            exit(EXIT_FAILURE);
        }
///STEP 5. Reset time to required restart time.
        _cController.initializerestart(_restart->getrestartime(),_minimizationTime);
	    #ifdef SLOWDOWNINITIALCYCLE
	    _slowedminimizationcutoffTime += _restart->getrestartime();
		#endif
        cout<<"Tau reset to "<<tau()<<endl;
///STEP 6. Reinitialize CBound eqlen, numHeads and numBoundHeads values as required by
/// datadump
        //sets
        _restart->CBoundinitializerestart();
///STEP 7. Assign copynumbers based on Chemistry input file or the datadump file as
// required by the user.
        _restart->restartupdateCopyNumbers();
        _restart->crosscheck();
        cout<<"Diffusion rates restored, diffusing molecules redistributed."<<endl;


//Step 8. re-add pin positions
        SystemParser p(_inputFile);
        FilamentSetup filSetup = p.readFilamentSetup();

        if(SysParams::Mechanics().pinBoundaryFilaments){
            PinRestartParser ppin(_inputDirectory + filSetup.pinRestartFile);
            ppin.resetPins();}

//Step 9. run mcontroller, update system, turn off restart state.
        updatePositions();
        updateNeighborLists();

        mins = chrono::high_resolution_clock::now();
        cout<<"Minimizing energy"<<endl;

        _subSystem.prevMinResult = _mController.run(false);
#ifdef OPTIMOUT
        mine= chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_runm(mine - mins);
        minimizationtime += elapsed_runm.count();
        std::cout<<"Time taken for minimization "<<elapsed_runm.count()<<endl;
#endif
        //DO NOT MOVE THIS LINE
        SysParams::RUNSTATE=true;

        //reupdate positions and neighbor lists
        updatePositions();
        updateNeighborLists();

//Step 10. Set Off rates back to original value.
        for(auto LL : Linker::getLinkers())
        {
            LL->getCLinker()->setOffRate(LL->getCLinker()->getOffReaction()->getBareRate());
            LL->updateReactionRates();
            LL->getCLinker()->getOffReaction()->updatePropensity();
            /*cout<<"L "<<LL->getId()<<" "<<LL->getMLinker()->getEqLength()<<" "
                << LL->getCLinker()->getOffRate()<<" "
                <<LL->getCLinker()->getOffReaction()->getBareRate()<<" "
                    <<LL->getMLinker()->stretchForce<<endl;*/

        }
        for(auto MM : MotorGhost::getMotorGhosts())
        {
            MM->getCMotorGhost()->setOffRate(MM->getCMotorGhost()->getOffReaction()->getBareRate());
            MM->updateReactionRates();
            MM->getCMotorGhost()->getOffReaction()->updatePropensity();
            /*cout<<"M "<<MM->getId()<<" "<<MM->getMMotorGhost()->getEqLength()<<" "
                << MM->getCMotorGhost()->getOffRate()<<" "
                <<MM->getCMotorGhost()->getOffReaction()->getBareRate()<<" "
                <<MM->getMMotorGhost()->stretchForce<<endl;*/
        }
        int dummy=0;
        for (auto BB: BranchingPoint::getBranchingPoints()) {
            dummy++;
            BB->getCBranchingPoint()->setOffRate(BB->getCBranchingPoint()->getOffReaction()->getBareRate());
            BB->updateReactionRates();
            BB->getCBranchingPoint()->getOffReaction()->updatePropensity();
            /*cout<<"B "<<BB->getId()<<" "<<BB->getMBranchingPoint()->getEqLength()<<" "
                << BB->getCBranchingPoint()->getOffRate()<<" "
                <<BB->getCBranchingPoint()->getOffReaction()->getBareRate()<<endl;*/
        }
//STEP 11: Get cylinders, activate filament reactions.
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                x->getCCylinder()->activatefilreactions();
                x->getCCylinder()->activatefilcrossreactions();
            }}

//Step 11b. Activate general reactions.
        for(auto C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto& rxn : C->getInternalReactionContainer().reactions()) {
                if(rxn->getReactionType() == ReactionType::REGULAR)
                    rxn->activateReaction();
            }}
        cout<<"Unbinding rates of bound species restored. filament reactions activated"<<endl;
//@
#ifdef CHEMISTRY
        _subSystem.updateBindingManagers();
#endif
#ifdef DYNAMICRATES
        updateReactionRates();
#endif
        delete _restart;
        cout<< "Restart procedures completed. Starting original Medyan framework"<<endl;
        cout << "---" << endl;
        cout << "Current simulation time = "<< tau() << endl;
        cout << endl;
        //restart phase ends

        //Crosscheck tau to make sure heap is ordered accurately.
        _cController.crosschecktau();
    }
#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    tauDatadump = tau();
    oldTau = 0;
#endif

#ifdef MECHANICS
    cout<<"Minimizing energy"<<endl;
    mins = chrono::high_resolution_clock::now();
    // update neighorLists before and after minimization. Need excluded volume
    // interactions.
	_subSystem.resetNeighborLists();
    auto minimizationResult = _mController.run(false);
    _subSystem.prevMinResult = minimizationResult;
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runm2(mine - mins);
    minimizationtime += elapsed_runm2.count();
#ifdef OPTIMOUT
    std::cout<<"Time taken for minimization "<<elapsed_runm2.count()<<endl;
#endif

    //activate/deactivate compartments
    mins = chrono::high_resolution_clock::now();
    //set initial values of variables.
    int tsaxis = SysParams::Boundaries().transfershareaxis;
    floatingpoint systemspan = 0.0;
    if(tsaxis == 0)
        systemspan = SysParams::Geometry().NX * SysParams::Geometry().compartmentSizeX;
    else if(tsaxis == 1)
        systemspan = SysParams::Geometry().NY * SysParams::Geometry().compartmentSizeY;
    else if(tsaxis == 2)
        systemspan = SysParams::Geometry().NZ * SysParams::Geometry().compartmentSizeZ;
    //copy vector to prevcopy
    bounds_prev[1] = systemspan;bounds_prev[0] = 0.0;
    bounds[1] = systemspan; bounds[0] =  0.0;
    activatedeactivateComp();
    moveBoundary(0.0);
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
    specialtime += elapsed_runspl.count();

    //reupdate positions and neighbor lists
    mins = chrono::high_resolution_clock::now();
    updatePositions();
#ifdef OPTIMOUT
    cout<<"Positions updated"<<endl;
#endif
    updateNeighborLists();
#ifdef OPTIMOUT
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runnl(mine - mins);
    nltime += elapsed_runnl.count();
    std::cout<<"NL time "<<elapsed_runnl.count()<<endl;
    mins = chrono::high_resolution_clock::now();
#endif
#ifdef DYNAMICRATES
    updateReactionRates();
#endif
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runrxn(mine - mins);
    rxnratetime += elapsed_runrxn.count();
#endif

#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    tauDatadump = tau();
    oldTau = 0;
#endif
    for(auto& o: _outputs) o->print(0);
    for(auto& o: _outputdump) o->print(0);

    resetCounters();

    cout << "Starting simulation..." << endl;

    int i = 1;

    //if runtime was specified, use this
    if(!areEqual(_runTime, 0.0)) {

#ifdef CHEMISTRY
        //activate/deactivate compartments
        mins = chrono::high_resolution_clock::now();
        activatedeactivateComp();
	    // set initial mechanical energy of system through a call to force field manager if dissipation tracking is enabled
	    if(SysParams::CParams.dissTracking){
		    _dt->setG1(minimizationResult.energiesAfter);
	    }
        mine= chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
        specialtime += elapsed_runspl.count();
        while(tau() <= _runTime) {
            auto minwhile = chrono::high_resolution_clock::now();
            //run ccontroller
            #ifdef OPTIMOUT
            cout<<"Starting chemistry"<<endl;
			#endif
            SysParams::DURINGCHEMISTRY = true;
            mins = chrono::high_resolution_clock::now();
            float factor = 1.0;
#ifdef SLOWDOWNINITIALCYCLE
            if(tau() <=_slowedminimizationcutoffTime)
            	factor = 10.0;
#endif
            floatingpoint chemistryTime = _minimizationTime/factor;
            #ifdef CROSSCHECK_CYLINDER
            string crosscheckchemname = _outputDirectory + "crosscheckChem.traj";
            if(CController::_crosscheckdumpFilechem.is_open())
                CController::_crosscheckdumpFilechem.close();
            CController::_crosscheckdumpFilechem.open(crosscheckchemname);
            #endif
            auto var = !_cController.run(chemistryTime);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runchem(mine - mins);
            chemistrytime += elapsed_runchem.count();
            SysParams::DURINGCHEMISTRY = false;
#ifdef OPTIMOUT
	        auto mtimex = CUDAcommon::tmin;
	        cout<<"motorbinding calls "<<mtimex.motorbindingcalls<<endl;
	        cout<<"motorunbinding calls "<<mtimex.motorunbindingcalls<<endl;
	        cout<<"motorwalking calls "<<mtimex.motorwalkingcalls<<endl;
	        cout<<"linkerbinding calls "<<mtimex.linkerbindingcalls<<endl;
	        cout<<"linkerunbinding calls "<<mtimex.linkerunbindingcalls<<endl;
	        CUDAcommon::tmin.motorbindingcalls = 0;
	        CUDAcommon::tmin.motorunbindingcalls = 0;
	        CUDAcommon::tmin.motorwalkingcalls = 0;
	        CUDAcommon::tmin.linkerbindingcalls = 0;
	        CUDAcommon::tmin.linkerunbindingcalls = 0;
#endif
            //print output if chemistry fails.
            mins = chrono::high_resolution_clock::now();
            if(var) {
                short counter = 0;
                for(auto& o: _outputs) { o->print(i); }
                resetCounters();
                break;
            }

            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runout(mine - mins);
            outputtime += elapsed_runout.count();
            //add the last step
            tauLastSnapshot += tau() - oldTau;
            tauLastMinimization += tau() - oldTau;
            tauLastNeighborList += tau() - oldTau;
            tauDatadump += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)

            //run mcontroller, update system
            if(tauLastMinimization >= _minimizationTime/factor) {
	            //check before rearrange
#ifdef CROSSCHECK_CYLINDER
                string crosscheckname = _outputDirectory + "crosscheckcyl.traj";
                Cylinder::_crosscheckdumpFile.open(crosscheckname);
                if(!Cylinder::_crosscheckdumpFile.is_open()) {
                    Cylinder::_crosscheckdumpFile << "There was an error opening file " << crosscheckname
                         << " for output. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                Cylinder::_crosscheckdumpFile << "Opening file " << crosscheckname << endl;
                Cylinder::_crosscheckdumpFile << "NCylinders " << Cylinder::getCylinders
                ().size() << endl;
#endif
                mins = chrono::high_resolution_clock::now();
                Bead::rearrange();
                Cylinder::updateAllData();

                string crosscheckmechname = _outputDirectory + "crosscheckmech.traj";
                CGMethod::_crosscheckdumpMechFile.open(crosscheckmechname);

                minimizationResult = _mController.run();
                _subSystem.prevMinResult = minimizationResult;
#ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile.close();
#endif
                mine= chrono::high_resolution_clock::now();

                
                chrono::duration<floatingpoint> elapsed_runm3(mine - mins);
                minimizationtime += elapsed_runm3.count();
                #ifdef OPTIMOUT
                std::cout<<"Time taken for minimization "<<elapsed_runm3.count()<<endl;
				#endif

                //update position
                mins = chrono::high_resolution_clock::now();
                updatePositions();
#ifdef CROSSCHECK_CYLINDER
                Cylinder::_crosscheckdumpFile.close();
#endif

                #ifdef OPTIMOUT
                cout<<"Position updated"<<endl;
				#endif

                tauLastMinimization = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_rxn2(mine - mins);
                updateposition += elapsed_rxn2.count();

                // perform multiple functions to update cumulative energy counters and reset the mechanical energy variables
                if(SysParams::CParams.dissTracking){
                    _dt->setGMid(minimizationResult.energiesBefore);
                    _dt->updateAfterMinimization(minimizationResult.energiesAfter);
                }

	            //update reaction rates
	            mins = chrono::high_resolution_clock::now();
#ifdef DYNAMICRATES
	            updateReactionRates();
#ifdef OPTIMOUT
	            cout<<"updated Reaction Rates"<<endl;
#endif
#endif
	            mine= chrono::high_resolution_clock::now();
	            chrono::duration<floatingpoint> elapsed_rxn3(mine - mins);
	            rxnratetime += elapsed_rxn3.count();
            }


            //output snapshot
            if(tauLastSnapshot >= _snapshotTime) {
                mins = chrono::high_resolution_clock::now();
                cout << "Current simulation time = "<< tau() << endl;
                for(auto& o: _outputs) o->print(i);
                resetCounters();
                i++;
                tauLastSnapshot = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_runout2(mine - mins);
                outputtime += elapsed_runout2.count();
#ifdef OPTIMOUT
                cout<< "Chemistry time for cycle=" << chemistrytime <<endl;
                cout << "Minimization time for cycle=" << minimizationtime <<endl;
                cout<< "Neighbor-list+Bmgr-time for cycle="<<nltime<<endl;
                cout<<"update-position time for cycle="<<updateposition<<endl;

                cout<<"rxnrate time for cycle="<<rxnratetime<<endl;
                cout<<"Output time for cycle="<<outputtime<<endl;
                cout<<"Special time for cycle="<<specialtime<<endl;
#endif
            }
            if(tauDatadump >= _datadumpTime) {
                for (auto& o: _outputdump) o->print(0);
                tauDatadump = 0.0;
            }
#elif defined(MECHANICS)
            for(auto& o: _outputs) o->print(i);
	        resetCounters();
            i++;
#endif

#ifdef CHEMISTRY
            //activate/deactivate compartments
            mins = chrono::high_resolution_clock::now();
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
            specialtime += elapsed_runspl.count();

            // update neighbor lists & Binding Managers
            if(tauLastNeighborList >= _neighborListTime/factor) {
                mins = chrono::high_resolution_clock::now();
                updateNeighborLists();
                tauLastNeighborList = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_runnl2(mine - mins);
                nltime += elapsed_runnl2.count();
#ifdef OPTIMOUT
                cout<<"update NeighborLists"<<endl;
#endif
            }
            //Special protocols
            mins = chrono::high_resolution_clock::now();
            //special protocols
            executeSpecialProtocols();
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runspl2(mine - mins);
            specialtime += elapsed_runspl2.count();
            oldTau = tau();
#ifdef CUDAACCL

            //reset CUDA context
//            CUDAcommon::handleerror(cudaDeviceSynchronize(), "cudaDeviceSynchronize", "Controller.cu");

//            CUDAcommon::handleerror(cudaDeviceReset(), "cudaDeviceReset", "Controller.cu");


//            size_t free, total;
//            CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
//            fprintf(stdout,"\t### After Reset Available VRAM : %g Mo/ %g Mo(total)\n\n",
//                    free/1e6, total/1e6);
//
//            cudaFree(0);
//
//            CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
//            fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
//                    free/1e6, total/1e6);
#ifdef CUDAACCL
                //@{
    size_t free2;
    CUDAcommon::handleerror(cudaMemGetInfo(&free2, &total));
    cudaFree(0);
    CUDAcommon::handleerror(cudaMemGetInfo(&free2, &total));
    std::cout<<"Free VRAM after CUDA operations in bytes "<<free2<<". Total VRAM in bytes "
             <<total<<endl;
    std::cout<<"Lost VRAM in bytes "<<CUDAcommon::getCUDAvars().memincuda-free2<<endl;
            std::cout<<endl;
    //@}
#endif
#endif
        }
#endif
    }
    //if run steps were specified, use this
    if(_runSteps != 0) {

#ifdef CHEMISTRY
        while(totalSteps <= _runSteps) {
            //run ccontroller
            if(!_cController.runSteps(_minimizationSteps)) {
                for(auto& o: _outputs) o->print(i);
                resetCounters();
                break;
            }

            //add the last step
            stepsLastSnapshot += _minimizationSteps;
            stepsLastMinimization += _minimizationSteps;
            stepsLastNeighborList += _minimizationSteps;

            totalSteps += _minimizationSteps;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            if(stepsLastMinimization >= _minimizationSteps) {
                Bead::rearrange();
                Cylinder::updateAllData();
                _mController.run();
                updatePositions();

#ifdef DYNAMICRATES
                updateReactionRates();
#endif

                stepsLastMinimization = 0;
            }

            if(stepsLastSnapshot >= _snapshotSteps) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto& o: _outputs) o->print(i);
                resetCounters();
                i++;
                stepsLastSnapshot = 0;
            }
#elif defined(MECHANICS)
            for(auto& o: _outputs) o->print(i);
            resetCounters();
            i++;
#endif

#ifdef CHEMISTRY
            //activate/deactivate compartments
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);

            // update neighbor lists
            if(stepsLastNeighborList >= _neighborListSteps) {
                updateNeighborLists();
                stepsLastNeighborList = 0;
            }

            //special protocols
            executeSpecialProtocols();
        }
#endif
    }

    //print last snapshots
    for(auto& o: _outputs) o->print(i);
	resetCounters();
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(chk2-chk1);
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
	#ifdef OPTIMOUT
    cout<< "Chemistry time for run=" << chemistrytime <<endl;
    cout << "Minimization time for run=" << minimizationtime <<endl;
    cout<< "Neighbor-list+Bmgr-time for run="<<nltime<<endl;
    cout<< "Neighbor-list time for run="<<nl2time<<endl;
    cout<< "Bmgr-vec time for run="<<bmgrvectime<<endl;
    cout<< "SIMD time for run="<<SubSystem::SIMDtime<<endl;
    cout<< "HYBD time for run="<<SubSystem::HYBDtime<<endl;
    cout<< "Bmgr time for run="<<bmgrtime<<endl;
    cout<<"update-position time for run="<<updateposition<<endl;
    cout<<"rxnrate time for run="<<rxnratetime<<endl;
    cout<<"Output time for run="<<outputtime<<endl;
    cout<<"Special time for run="<<specialtime<<endl;
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout<<"-----------"<<endl;
    if(true){
        auto mtime = CUDAcommon::tmin;
        cout<<"update-position time for run="<<updateposition<<endl;
        cout<<"update-position-cylinder time for run="<<updatepositioncylinder<<endl;
        cout<<"update-position-movable time for run="<<updatepositionmovable<<endl;
        cout<<"move-compartment cylinder ="<<mtime.timecylinderupdate<<" calls "<<mtime
        .callscylinderupdate<<endl;
        cout<<"move-compartment linker ="<<mtime.timelinkerupdate<<" calls "<<mtime
                .callslinkerupdate<<endl;
        cout<<"move-compartment motor ="<<mtime.timemotorupdate<<" calls "<<mtime
                .callsmotorupdate<<endl;
        cout<<"-----------"<<endl;

        cout << "Minimization time for run=" << minimizationtime <<endl;
        cout<<"Printing minimization times in seconds."<<endl;
        cout<<"starting minimization "<<mtime.vectorize<<endl;
        cout<<"Finding lambda "<<mtime.findlambda<<endl;
        cout<<"copy forces "<<mtime.copyforces<<endl;
        cout<<"other computations "<<mtime.tother<<endl;
        cout<<"end minimization "<<mtime.endminimization<<endl;
        cout<<"compute energies "<<mtime.computeenergy<<" calls "<<mtime
        .computeenergycalls<<endl;
        cout<<"compute energieszero "<<mtime.computeenergyzero<<" calls "<<mtime
        .computeenerycallszero<<endl;
	    cout<<"compute energiesnonzero "<<mtime.computeenergynonzero<<" calls "<<mtime
			    .computeenerycallsnonzero<<endl;
        cout<<"compute forces "<<mtime.computeforces<<" calls "<<mtime
                .computeforcescalls<<endl;
        cout<<"Time taken to compute energy in each forcefield "<<endl;
        cout<<"Filament, Linker, Motor, Branching, Excluded Volume, and Boundary"<<endl;
        for(auto x:mtime.individualenergies)
            cout<<x<<" ";
        cout<<endl;
	    cout<<"Time taken to compute energy in "<<endl;
	    cout<<"Filament Stretching "<<mtime.stretchingenergy<<endl;
	    cout<<"Filament Bending "<<mtime.bendingenergy<<endl;
        cout<<"Time taken to compute energyzero in each forcefield "<<endl;
        for(auto x:mtime.individualenergieszero)
            cout<<x<<" ";
        cout<<endl;
        cout<<"Time taken to compute energynonzero in each forcefield "<<endl;
        for(auto x:mtime.individualenergiesnonzero)
            cout<<x<<" ";
        cout<<endl;
        cout<<"Time taken to compute forces in each forcefield "<<endl;
        for(auto x:mtime.individualforces)
            cout<<x<<" ";
        cout<<endl;
	    cout<<"Time taken to compute forces in "<<endl;
	    cout<<"Filament Stretching "<<mtime.stretchingforces<<endl;
	    cout<<"Filament Bending "<<mtime.bendingforces<<endl;

		cout<<"Number of interactions considered in each force field"<<endl;

		cout<<"Filament Stretching "<<mtime.numinteractions[0]<<endl;
	    cout<<"Filament Bending "<<mtime.numinteractions[1]<<endl;
	    cout<<"Linker Stretching "<<mtime.numinteractions[2]<<endl;
	    cout<<"Motor Stretching "<<mtime.numinteractions[3]<<endl;
	    cout<<"Branching Stretching "<<mtime.numinteractions[4]<<endl;
	    cout<<"Branching Bending "<<mtime.numinteractions[5]<<endl;
	    cout<<"Branching Dihedral "<<mtime.numinteractions[6]<<endl;
	    cout<<"Branching Position "<<mtime.numinteractions[7]<<endl;
	    cout<<"Cylinder-Cylinder Repulsion "<<mtime.numinteractions[8]<<endl;
	    cout<<"Cylinder-Boundary Repulsion "<<mtime.numinteractions[9]<<endl;
    }
    if(true) {
        cout << "Printing callback times" << endl;
        auto ctime = CUDAcommon::ctime;
        auto ccount = CUDAcommon::ccount;
        cout << "UpdateBrancherBindingCallback " << ctime.tUpdateBrancherBindingCallback
             << " count "
             << ccount.cUpdateBrancherBindingCallback << endl;
        cout << "UpdateLinkerBindingCallback " << ctime.tUpdateLinkerBindingCallback
             << " count "
             << ccount.cUpdateLinkerBindingCallback << endl;
        cout << "UpdateMotorBindingCallback " << ctime.tUpdateMotorBindingCallback
             << " count "
             << ccount.cUpdateMotorBindingCallback << endl;
        cout << "UpdateMotorIDCallback " << ctime.tUpdateMotorIDCallback << " count "
             << ccount.cUpdateMotorIDCallback << endl;
        cout << "FilamentExtensionPlusEndCallback "
             << ctime.tFilamentExtensionPlusEndCallback << " count "
             << ccount.cFilamentExtensionPlusEndCallback << endl;
        cout << "FilamentExtensionMinusEndCallback "
             << ctime.tFilamentExtensionMinusEndCallback << " count "
             << ccount.cFilamentExtensionMinusEndCallback << endl;
        cout << "FilamentRetractionPlusEndCallback "
             << ctime.tFilamentRetractionPlusEndCallback << " count "
             << ccount.cFilamentRetractionPlusEndCallback << endl;
        cout << "FilamentRetractionMinusEndCallback "
             << ctime.tFilamentRetractionMinusEndCallback << " count "
             << ccount.cFilamentRetractionMinusEndCallback << endl;
        cout << "FilamentPolymerizationPlusEndCallback "
             << ctime.tFilamentPolymerizationPlusEndCallback << " count "
             << ccount.cFilamentPolymerizationPlusEndCallback << endl;
        cout << "FilamentPolymerizationMinusEndCallback "
             << ctime.tFilamentPolymerizationMinusEndCallback << " count "
             << ccount.cFilamentPolymerizationMinusEndCallback << endl;
        cout << "FilamentDepolymerizationPlusEndCallback "
             << ctime.tFilamentDepolymerizationPlusEndCallback << " count "
             << ccount.cFilamentDepolymerizationPlusEndCallback << endl;
        cout << "FilamentDepolymerizationMinusEndCallback "
             << ctime.tFilamentDepolymerizationMinusEndCallback << " count "
             << ccount.cFilamentDepolymerizationMinusEndCallback << endl;
        cout << "BranchingPointUnbindingCallback " << ctime.tBranchingPointUnbindingCallback
             << " count "
             << ccount.cBranchingPointUnbindingCallback << endl;
        cout << "BranchingCallback " << ctime.tBranchingCallback << " count "
             << ccount.cBranchingCallback << endl;
        cout << "LinkerUnbindingCallback " << ctime.tLinkerUnbindingCallback << " count "
             << ccount.cLinkerUnbindingCallback << endl;
        cout << "LinkerBindingCallback " << ctime.tLinkerBindingCallback << " count "
             << ccount.cLinkerBindingCallback << endl;
        cout << "MotorUnbindingCallback " << ctime.tMotorUnbindingCallback << " count "
             << ccount.cMotorUnbindingCallback << endl;
        cout << "MotorBindingCallback " << ctime.tMotorBindingCallback << " count "
             << ccount.cMotorBindingCallback << endl;
        cout << "MotorWalkingCallback " << ctime.tMotorWalkingCallback << " count "
             << ccount.cMotorWalkingCallback << endl;
        cout << "MotorMovingCylinderCallback " << ctime.tMotorMovingCylinderCallback
             << " count "
             << ccount.cMotorMovingCylinderCallback << endl;
        cout << "FilamentCreationCallback " << ctime.tFilamentCreationCallback << " count "
             << ccount.cFilamentCreationCallback << endl;
        cout << "FilamentSeveringCallback " << ctime.tFilamentSeveringCallback << " count "
             << ccount.cFilamentSeveringCallback << endl;
        cout << "FilamentDestructionCallback " << ctime.tFilamentDestructionCallback
             << " count "
             << ccount.cFilamentDestructionCallback << endl;
        cout << "------------" << endl;
        cout << "Printing neighbor times" << endl;
        cout << "Dynamic neighbor " << SubSystem::timedneighbor << endl;
        cout << "Neighbor " << SubSystem::timeneighbor << endl;
        cout << "Trackable " << SubSystem::timetrackable << endl;

        cout << "-------------" << endl;
        cout << "Filament extendPlusEnd 1 " << Filament::FilextendPlusendtimer1 << endl;
        cout << "Filament extendPlusEnd 2 " << Filament::FilextendPlusendtimer2 << endl;
        cout << "-------------" << endl;
        cout << "Cylinder constructor" << endl;
        cout << "part1 " << Cylinder::timecylinder1 << " part2 " << Cylinder::timecylinder2
             << " "
                "Ccylinder "
             << Cylinder::timecylinderchem << " mCylinder " << Cylinder::timecylindermech
             << endl;
        cout << "initializeCCylinder for loop " << ChemManager::tchemmanager1 << endl;
        cout << "extension Front/Back " << ChemManager::tchemmanager2 << endl;
        cout << "initialize " << ChemManager::tchemmanager3 << endl;
        cout << "last part " << ChemManager::tchemmanager4 << endl;
        cout << "------------" << endl;
        cout << "PolyPlusEndTemplate time" << endl;
        cout << "For loop " << CUDAcommon::ppendtime.rxntempate1 << " part2 (findspecies) "
             << CUDAcommon::ppendtime.rxntempate2 << " part3 (create rxn) "
             << CUDAcommon::ppendtime
                     .rxntempate3 << " part4 (Callback) "
             << CUDAcommon::ppendtime.rxntempate4 << endl;
        cout<<" Displaying chemistry times"<<endl;
        cout<<"Counts fired for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.reactioncount[i]<<" ";
        cout<<endl;
        cout<<"Time taken to fire each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.totaltime[i]<<" ";
        cout<<endl;
        cout<<"Time taken to emitSignal for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.emitsignal[i]<<" ";
        cout<<endl;
        cout<<"Time taken for dependency updates for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.dependencytime[i]<<" ";
        cout<<endl;
        cout<<"Total number of dependencies for reactions fired based on ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.dependentrxncount[i]<<" ";
        cout<<endl;
    }
	#endif
    cout << "Done with simulation!" << endl;
#ifdef CUDAACCL
    cudaDeviceReset();
#endif
}
