
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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
#include <unordered_set>

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
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneHierarchy.h"
#include "Structure/SurfaceMesh/MembraneRegion.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "MController.h"
#include "Cylinder.h"
#include <unordered_map>
#include <tuple>
#include <vector>
#include <algorithm>
#include "Rand.h"
#include "Restart.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "util/io/log.h"
#include "util/profiler.h"
#include "Visual/Window.hpp"
using namespace mathfunc;

namespace {

void rearrangeAllDatabases() {
    Cylinder::rearrange(); Cylinder::updateData();
}

} // namespace

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
    _outputs.push_back(new BasicSnapshot(_outputDirectory + "snapshot.traj", _subSystem));
    _outputs.push_back(new BirthTimes(_outputDirectory + "birthtimes.traj", _subSystem));
    _outputs.push_back(new Forces(_outputDirectory + "forces.traj", _subSystem));
    _outputs.push_back(new Tensions(_outputDirectory + "tensions.traj", _subSystem));

    _outputs.push_back(new PlusEnd(_outputDirectory + "plusend.traj", _subSystem));
    //ReactionOut should be the last one in the output list
    //Otherwise incorrect deltaMinusEnd or deltaPlusEnd values may be genetrated.
    _outputs.push_back(new ReactionOut(_outputDirectory + "monomers.traj", _subSystem));
    //add br force out and local diffussing species concentration
    //_outputs.push_back(new BRForces(_outputDirectory + "repulsion.traj", _subSystem));
    //_outputs.push_back(new PinForces(_outputDirectory + "pinforce.traj", _subSystem));

    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters()) exit(EXIT_FAILURE);

    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "---" << endl;
    LOG(STEP) << "Initializing geometry...";
    _gController->initializeGrid();
    LOG(INFO) << "Done.";
    
    //Initialize boundary
    cout << "---" << endl;
    LOG(STEP) << "Initializing boundary...";
    
    auto BTypes = p.readBoundaryType();
    p.readBoundParams();

    //initialize
    _gController->initializeBoundary(BTypes);
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
    _mController->initialize(MTypes, MAlgorithm);
    LOG(INFO) << "Done.";

#endif

#ifdef CHEMISTRY
    //Activate necessary compartments for diffusion
    _gController->setActiveCompartments();

    if(_subSystem->getBoundary()->getShape() == BoundaryShape::Cylinder){
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
            C->computeSlicedVolumeArea(Compartment::SliceMethod::CylinderBoundary);
        }
    }
    else{
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
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
    _cController->initialize(CAlgorithm.algorithm, ChemData);
    LOG(INFO) << "Done.";
    
    //Set up chemistry output if any
    string chemsnapname = _outputDirectory + "chemistry.traj";
    _outputs.push_back(new Chemistry(chemsnapname, _subSystem, ChemData,
                                     _subSystem->getCompartmentGrid()));

    string concenname = _outputDirectory + "concentration.traj";
    _outputs.push_back(new Concentrations(concenname, _subSystem, ChemData));

#endif

#ifdef DYNAMICRATES
    cout << "---" << endl;
    LOG(STEP) << "Initializing dynamic rates...";
    //read dynamic rate parameters
    p.readDyRateParams();

    //read dynamic rate types
    DynamicRateType DRTypes = p.readDynamicRateType();

    //init controller
    _drController->initialize(DRTypes);
    LOG(INFO) << "Done.";
    
#endif

    //Check consistency of all chemistry and mechanics parameters
    cout << "---" << endl;
    LOG(STEP) << "Checking cross-parameter consistency...";
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

    // Initialize the membrane mesh adapter
    // Currently the values are simply represented as magic numbers
    _meshAdapter = std::make_unique< adaptive_mesh::MembraneMeshAdapter >(
        typename adaptive_mesh::MembraneMeshAdapter::Parameter {
            // Topology
            adaptive_mesh::surface_mesh_min_degree,
            adaptive_mesh::surface_mesh_max_degree,
            adaptive_mesh::edge_flip_min_dot_normal,
            adaptive_mesh::edge_collapse_min_quality_improvement,
            adaptive_mesh::edge_collapse_min_dot_normal,
            // Relocation
            adaptive_mesh::vertex_relaxation_epsilon,
            adaptive_mesh::vertex_relaxation_dt,
            adaptive_mesh::vertex_relocation_max_iter_relocation,
            adaptive_mesh::vertex_relocation_max_iter_tot,
            // Size diffusion
            adaptive_mesh::size_measure_curvature_resolution,
            adaptive_mesh::size_measure_max,
            adaptive_mesh::size_measure_diffuse_iter,
            // Main loop
            adaptive_mesh::mesh_adaptation_topology_max_iter,
            adaptive_mesh::mesh_adaptation_soft_max_iter,
            adaptive_mesh::mesh_adaptation_hard_max_iter
        }
    );
    
    LOG(INFO) << "Done.";
    
    //setup initial network configuration
    setupInitialNetwork(p);
#ifdef HYBRID_NLSTENCILLIST
//    _subSystem->initializeHNeighborList();
#endif
    //setup special structures
    setupSpecialStructures(p);

    SysParams::INITIALIZEDSTATUS = true;
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

    /**************************************************************************
    Now starting to add the membrane into the network.
    **************************************************************************/
    MembraneSetup MemSetup = p.readMembraneSetup();
    
    cout << "---" << endl;
    cout << "Initializing membranes...";

    std::vector< MembraneParser::MembraneInfo > membraneData;
    if(MemSetup.inputFile != "") {
        membraneData = MembraneParser(_inputDirectory + MemSetup.inputFile).readMembranes();
    }
    // Membrane auto initializer is currently not provided,
    // which means that the input file is the only source of membrane information.
    
    // add membranes
    for (auto& it: membraneData) {
        
        short type = 0; // Currently set as default(0).
        
        if(type >= SysParams::Chemistry().numMembranes) {
            cout << "Membrane data specified contains an invalid membrane type. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        Membrane* newMembrane = _subSystem->addTrackable<Membrane>(
            _subSystem,
            type,
            it.vertexCoordinateList,
            it.triangleVertexIndexList
        );

        // Update membrane geometry
        newMembrane->updateGeometryValue();

        // Add to the global membrane hierarchy
        MembraneHierarchy::addMembrane(newMembrane);
    }
    cout << "Done. " << membraneData.size() << " membranes created." << endl;

    // Create a region inside the membrane
    auto regionInMembrane = (
        membraneData.size()?
        MembraneRegion::makeByChildren(MembraneHierarchy::getRoot()):
        make_unique<MembraneRegion>(_subSystem->getBoundary())
    );

    // Optimize the membrane
    membraneAdaptiveRemesh();

    // Deactivate all the compartments outside membrane, and mark boundaries as interesting
    for(auto c : _subSystem->getCompartmentGrid()->getCompartments()) {
        if(!c->getTriangles().empty()) {
            // Contains triangles, so this compartment is at the boundary.
            c->boundaryInteresting = true;

            // Update partial activate status
            c->computeSlicedVolumeArea(Compartment::SliceMethod::Membrane);
            _cController->updateActivation(c, Compartment::ActivateReason::Membrane);

        } else if( ! regionInMembrane->contains(vector2Vec<3, double>(c->coordinates()))) {
            // Compartment is outside the membrane
            _cController->deactivate(c, true);
        }
    }

    // Transfer species from all the inactive compartments
    {
        vector<Compartment*> ac, ic;
        for(auto c : _subSystem->getCompartmentGrid()->getCompartments()) {
            if(c->isActivated()) ac.push_back(c);
            else                 ic.push_back(c);
        }
        auto nac = ac.size();
        for(auto c : ic) {
            for(auto &sp : c->getSpeciesContainer().species()) {
                int copyNumber = sp->getN();
                unordered_set<Species*> sp_targets;
                if(sp->getFullName().find("Bound") == string::npos) {
                    while(copyNumber > 0) {
                        sp->down();
                        auto tc = ac[Rand::randInteger(0, nac-1)];
                        auto sp_target = tc->findSpeciesByName(sp->getName());
                        sp_targets.insert(sp_target);
                        sp_target->up();
                        --copyNumber;
                    }
                }
                for(auto sp_target : sp_targets)
                    sp_target->updateReactantPropensities();
                sp->updateReactantPropensities();
            }
        }
    }

    /**************************************************************************
    Now starting to add the filaments into the network.
    **************************************************************************/
    // Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
    
    cout << "---" << endl;
    cout << "Initializing filaments...";

    if(FSetup.inputFile != "") {
        FilamentParser fp(_inputDirectory + FSetup.inputFile);
        filaments = fp.readFilaments();
    }
    fil=get<0>(filaments);
    //add other filaments if specified
    FilamentInitializer* fInit = new RandomFilamentDist();
    
    auto filamentsGen = fInit->createFilaments(*regionInMembrane,
                                               FSetup.numFilaments,
                                               FSetup.filamentType,
                                               FSetup.filamentLength);
    auto filGen=get<0>(filamentsGen);
    fil.insert(fil.end(), filGen.begin(), filGen.end());
    delete fInit;

    //add filaments
    for (auto it: fil) {

        auto coord1 = get<1>(it);
        auto coord2 = get<2>(it);
        auto type = get<0>(it);

        if(type >= SysParams::Chemistry().numFilaments) {
            cout << "Filament data specified contains an "
                 <<"invalid filament type. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        vector<vector<double>> coords = {coord1, coord2};

        if(coord2.size()==3){

            double d = twoPointDistance(coord1, coord2);
            vector<double> tau = twoPointDirection(coord1, coord2);
            int numSegment = static_cast<int>(std::round(d / SysParams::Geometry().cylinderSize[type]));

            // check how many segments can fit between end-to-end of the filament
            if (numSegment == 0)
                _subSystem->addTrackable<Filament>(_subSystem, type, coords, 2, FSetup.projectionType);
            else
                _subSystem->addTrackable<Filament>(_subSystem, type, coords, numSegment + 1, FSetup.projectionType);
        }
        else if(coord2.size()>3){
            int numSegment = coord2.size()/3;
            vector<vector<double>> coords;
            coords.push_back(coord1);
            for(int id=0;id<numSegment;id++)
                coords.push_back({coord2[id*3],coord2[id*3+1],coord2[id*3+2]});

            if (numSegment == 0)
                _subSystem->addTrackable<Filament>(_subSystem, type, coords, 2, FSetup.projectionType);
            else
                _subSystem->addTrackable<Filament>(_subSystem, type, coords, numSegment + 1, FSetup.projectionType);
        }
    }
    cout << "Done. " << fil.size() << " filaments created." << endl;
    cout<<"Total cylinders "<<Cylinder::getCylinders().size()<<endl;
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
        filamentData fil=get<0>(filaments);
        for (auto it: fil) {

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

void Controller::activatedeactivateComp(){

    if(SysParams::Boundaries().transfershareaxis!=3){
        fCompmap.clear();
        bCompmap.clear();
        activatecompartments.clear();
        ControlfrontbackEndComp();
//            std::cout<<fCompmap.size()<<" "<<bCompmap.size()<<" "<<activatecompartments.size()<<endl;
        for(auto it=activatecompartments.begin();it!=activatecompartments.end();it++)
        {
            if(!(*it)->isActivated())
                _cController->activate(*it);
        }
        //deactivate compartments starting from the right extreme
        for (std::multimap<int,Compartment*>::reverse_iterator it=fCompmap.rbegin(); it!=fCompmap.rend(); ++it)
            _cController->deactivate(it->second);
        //deactivate compartments starting from the left extreme
        for (std::multimap<int,Compartment*>::iterator it=bCompmap.begin(); it!=bCompmap.end(); ++it)
            _cController->deactivate(it->second);
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
    double systemspan = 0.0;
    double cmpsize = 0.0;
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
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
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
        for (auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
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
        for (auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
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

void Controller::moveBoundary(double deltaTau) {
    //calculate distance to move
    double dist = SysParams::Boundaries().moveSpeed * deltaTau;

    if(SysParams::Boundaries().transfershareaxis != 3){
        vector<double> distvec= {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
        _subSystem->getBoundary()->move(distvec);
    }
    //calculate system volume.
    _subSystem->getBoundary()->volume();
}

void Controller::updateActiveCompartments() {
    // For this function to work, we must assume that each minimization step
    // will push the membrane boundary no more than 1 compartment, so that
    // changes will only happen at the neighborhood of the previous boundary.
    auto& allMembranes = Membrane::getMembranes();

    // Currently only the 0th membrane will be considered
    if(allMembranes.size()) {
        Membrane* theMembrane = allMembranes[0];
        // For non empty compartments, we mark them as interesting and update their status
        // For the "interesting" compartments last round but now empty, we fully activate or deactivate them
        // For the rest we do nothing, assuming that the membranes will NOT move across a whole compartment
        for(auto c: GController::getCompartmentGrid()->getCompartments()) {
            auto& ts = c->getTriangles();
            if(!ts.empty()) {
                // Update partial activate status
                c->computeSlicedVolumeArea(Compartment::SliceMethod::Membrane);
                _cController->updateActivation(c, Compartment::ActivateReason::Membrane);

                // No matter whether the compartment is interesting before, mark it as interesting
                c->boundaryInteresting = true;
            } else if(c->boundaryInteresting) { // Interesting last round but now empty
                bool inMembrane = (
                    (!theMembrane->isClosed()) ||
                    (theMembrane->contains(vector2Vec<3, double>(c->coordinates())))
                );
                if(inMembrane) {
                    // Fully activate the compartment
                    c->resetVolumeFrac();
					const auto& fullArea = GController::getCompartmentArea();
                    c->setPartialArea({{
                        fullArea[0], fullArea[0],
						fullArea[1], fullArea[1],
						fullArea[2], fullArea[2]
                    }});
                    _cController->updateActivation(c, Compartment::ActivateReason::Membrane);
                } else {
                    // Deactivate the compartment
                    _cController->deactivate(c);
                }

                // Mark the compartment as not interesting
                c->boundaryInteresting = false;
            }
        }
    } // Otherwise, no membrane exists. Do nothing
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
    chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();
    //vectorize cylinder to have all cylinder information in a few arrays.
    mine = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runbvec(mine - mins);
    bmgrvectime += elapsed_runbvec.count();
    mins = chrono::high_resolution_clock::now();
    //Full reset of neighbor lists
    _subSystem->resetNeighborLists();
    mine = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runnl2(mine - mins);
    nl2time += elapsed_runnl2.count();
#ifdef CHEMISTRY
    mins = chrono::high_resolution_clock::now();
    _subSystem->updateBindingManagers();
    mine = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runb(mine - mins);
    bmgrtime += elapsed_runb.count();
//    std::cout<<"time split "<<elapsed_runnl2.count()<<" "<<elapsed_runbvec.count()<<" "
//            ""<<elapsed_runb.count()<<endl;
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

            cout << _subSystem->getBoundary()->distance(b->vcoordinate()) << endl;
            cout << SysParams::Mechanics().pinDistance << endl;


            //if within dist to boundary, add
            if(_subSystem->getBoundary()->distance(b->vcoordinate()) < SysParams::Mechanics().pinDistance) {

                b->pinnedPosition = b->vcoordinate();
                b->addAsPinned();
            }
        }
    }
}
//Qin
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

            auto index = Rand::randDouble(0,1);
            //cout << index <<endl;
            //if within dist to boundary and index > 0.5, add
            if(_subSystem->getBoundary()->lowerdistance(b->vcoordinate()) < SysParams::Mechanics().pinDistance
               && index < SysParams::Mechanics().pinFraction && b->isPinned() == false) {
                //cout << index << endl;
                b->pinnedPosition = b->vcoordinate();
                b->addAsPinned();
            }
        }
    }
}

void Controller::membraneAdaptiveRemesh() const {
    // Requires _meshAdapter to be already initialized
    for(auto m : Membrane::getMembranes()) {
        _meshAdapter->adapt(m->getMesh());

        // Update necessary geometry for the system
        m->updateGeometryValueForSystem();
    }
}

void Controller::run() {

#ifdef CHEMISTRY
    double tauLastSnapshot = 0;
    double tauLastMinimization = 0;
    double tauLastNeighborList = 0;
    double oldTau = 0;

    long stepsLastSnapshot = 0;
    long stepsLastMinimization = 0;
    long stepsLastNeighborList = 0;

    long totalSteps = 0;
#endif
    chrono::high_resolution_clock::time_point chk1, chk2, mins, mine;
    chk1 = chrono::high_resolution_clock::now();
//RESTART PHASE BEGINS
    if(SysParams::RUNSTATE==false){
        cout<<"RESTART PHASE BEINGS."<<endl;
        Restart* _restart = new Restart(_subSystem, filaments,_chemData);
//Step 1. Turn off diffusion, passivate filament reactions and empty binding managers.
//        _restart->settorestartphase();
        cout<<"Turned off Diffusion, filament reactions."<<endl;
//Step 2. Add bound species to their respective binding managers. Turn off unbinding, update propensities.
        //_restart->addtoHeaplinkermotor();
        _restart->addtoHeapbranchers();
        _restart->addtoHeaplinkermotor();
        cout<<"Bound species added to reaction heap."<<endl;
//Step 2A. Turn off diffusion, passivate filament reactions and empty binding managers.
        _restart->settorestartphase();
//Step 3. ############ RUN LINKER/MOTOR REACTIONS TO BIND BRANCHERS, LINKERS, MOTORS AT RESPECTIVE POSITIONS.#######
        cout<<"Reactions to be fired "<<_restart->getnumchemsteps()<<endl;
        _cController->runSteps(_restart->getnumchemsteps());
        cout<<"Reactions fired! Displaying heap"<<endl;
//Step 4. Display the number of reactions yet to be fired. Should be zero.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto &Mgr:C->getFilamentBindingManagers()){
                int numsites = 0;
#ifdef NLORIGINAL
                numsites = Mgr->numBindingSites();
#endif
#ifdef NLSTENCILLIST
                numsites = Mgr->numBindingSitesstencil();
#endif
                if(numsites == 0)
                    cout<< numsites<<" ";
                else{
                    cout<<endl;
                    cout<<"Few reactions are not fired! Cannot restart this trajectory. Exiting ..."<<endl;
                    exit(EXIT_FAILURE);
                }
            }}
        cout<<endl;
        _restart->redistributediffusingspecies();
        cout<<"Diffusion rates restored, diffusing molecules redistributed."<<endl;

//Step 4.5. re-add pin positions
        SystemParser p(_inputFile);
        FilamentSetup filSetup = p.readFilamentSetup();

        if(SysParams::Mechanics().pinBoundaryFilaments){
            PinRestartParser ppin(_inputDirectory + filSetup.pinRestartFile);
            ppin.resetPins();}

//Step 5. run mcontroller, update system, turn off restart state.
        updatePositions();
        updateNeighborLists();

        mins = chrono::high_resolution_clock::now();
        cout<<"Minimizing energy"<<endl;

        _mController->run(false);

        mine= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runm(mine - mins);
        minimizationtime += elapsed_runm.count();
        std::cout<<"Time taken for minimization "<<elapsed_runm.count()<<endl;
        SysParams::RUNSTATE=true;

        //reupdate positions and neighbor lists
        updatePositions();
        updateNeighborLists();

//Step 6. Set Off rates back to original value.
        for(auto LL : Linker::getLinkers())
        {
            LL->getCLinker()->setOffRate(LL->getCLinker()->getOffReaction()->getBareRate());
            LL->getCLinker()->getOffReaction()->setRate(LL->getCLinker()->getOffReaction()->getBareRate());
            LL->updateReactionRates();
            LL->getCLinker()->getOffReaction()->updatePropensity();

        }
        for(auto MM : MotorGhost::getMotorGhosts())
        {
            MM->getCMotorGhost()->setOffRate(MM->getCMotorGhost()->getOffReaction()->getBareRate());
            MM->getCMotorGhost()->getOffReaction()->setRate(MM->getCMotorGhost()->getOffReaction()->getBareRate());
            MM->updateReactionRates();
            MM->getCMotorGhost()->getOffReaction()->updatePropensity();
        }
        int dummy=0;
        for (auto BB: BranchingPoint::getBranchingPoints()) {
            dummy++;
            BB->getCBranchingPoint()->setOffRate(BB->getCBranchingPoint()->getOffReaction()->getBareRate());
            BB->getCBranchingPoint()->getOffReaction()->setRate(BB->getCBranchingPoint()->getOffReaction()->getBareRate());
            BB->getCBranchingPoint()->getOffReaction()->updatePropensity();
        }
//STEP 7: Get cylinders, activate filament reactions.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                x->getCCylinder()->activatefilreactions();
                x->getCCylinder()->activatefilcrossreactions();
            }}
        cout<<"Unbinding rates of bound species restored. filament reactions activated"<<endl;
//@
#ifdef CHEMISTRY
        _subSystem->updateBindingManagers();
#endif
#ifdef DYNAMICRATES
        updateReactionRates();
#endif
        cout<< "Restart procedures completed. Starting original Medyan framework"<<endl;
        cout << "---" << endl;
        resetglobaltime();
        _cController->restart();
        cout << "Current simulation time = "<< tau() << endl;
        //restart phase ends
    }

    //perform first minimization
//#ifdef MECHANICS
//    _mController->run(false);
//
//    //reupdate positions and neighbor lists
//    updatePositions();
//    updateNeighborLists();
//
//#ifdef DYNAMICRATES
//    updateReactionRates();
//#endif
//
//#endif

#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    oldTau = 0;
#endif

#ifdef MECHANICS
    cout<<"Minimizing energy"<<endl;
    mins = chrono::high_resolution_clock::now();
    _mController->run(false);
    membraneAdaptiveRemesh();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runm2(mine - mins);
    minimizationtime += elapsed_runm2.count();
    std::cout<<"Time taken for minimization "<<elapsed_runm2.count()<<endl;

    //activate/deactivate compartments
    mins = chrono::high_resolution_clock::now();
    //set initial values of variables.
    int tsaxis = SysParams::Boundaries().transfershareaxis;
    double systemspan = 0.0;
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
    chrono::duration<double> elapsed_runspl(mine - mins);
    specialtime += elapsed_runspl.count();

    //reupdate positions and neighbor lists
    mins = chrono::high_resolution_clock::now();
    updatePositions();
    rearrangeAllDatabases();
    updateNeighborLists();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runnl(mine - mins);
    nltime += elapsed_runnl.count();
    std::cout<<"NL time "<<elapsed_runnl.count()<<endl;
    mins = chrono::high_resolution_clock::now();
#ifdef DYNAMICRATES
    updateReactionRates();
#endif
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runrxn(mine - mins);
    rxnratetime += elapsed_runrxn.count();
#endif

#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    oldTau = 0;
#endif
    for(auto o: _outputs) o->print(0);
    resetCounters();

    cout << "Starting simulation..." << endl;


    int i = 1;

    //if runtime was specified, use this
    if(!areEqual(_runTime, 0.0)) {

#ifdef CHEMISTRY
        //activate/deactivate compartments
        mins = chrono::high_resolution_clock::now();
        //activatedeactivateComp();
        mine= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runspl(mine - mins);
        specialtime += elapsed_runspl.count();

        while(tau() <= _runTime) {
            //run ccontroller
            mins = chrono::high_resolution_clock::now();
            auto var = !_cController->run(_minimizationTime);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_runchem(mine - mins);
            chemistrytime += elapsed_runchem.count();

            //print output if chemistry fails.
            mins = chrono::high_resolution_clock::now();
            if(var) {
                for(auto o: _outputs) o->print(i);
                resetCounters();
                break;
            }
            mine= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_runout(mine - mins);
            outputtime += elapsed_runout.count();

            //add the last step
            tauLastSnapshot += tau() - oldTau;
            tauLastMinimization += tau() - oldTau;
            tauLastNeighborList += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
#ifdef CUDAACCL
            //@{
    size_t free, total;
    CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
    cudaFree(0);
    CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
    std::cout<<"Free VRAM before CUDA operations in bytes "<<free<<". Total VRAM in bytes "
             <<total<<endl;
            auto cvars = CUDAcommon::getCUDAvars();
            cvars.memincuda = free;
            CUDAcommon::cudavars = cvars;
    //@}
#endif
            //run mcontroller, update system
            if(tauLastMinimization >= _minimizationTime) {

                mins = chrono::high_resolution_clock::now();
                _mController->run();

                // Membrane remeshing
                membraneAdaptiveRemesh();

                mine= chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_runm3(mine - mins);
                minimizationtime += elapsed_runm3.count();
//                std::cout<<"Time taken for minimization "<<elapsed_runm3.count()<<endl;
                //update position
                mins = chrono::high_resolution_clock::now();

                updatePositions();

                // Update activation of the compartments
                updateActiveCompartments();

                tauLastMinimization = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_rxn2(mine - mins);
                rxnratetime += elapsed_rxn2.count();

            }
            //output snapshot
            if(tauLastSnapshot >= _snapshotTime) {
                mins = chrono::high_resolution_clock::now();
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                resetCounters();
                i++;
                tauLastSnapshot = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_runout2(mine - mins);
                outputtime += elapsed_runout2.count();
            }

#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            resetCounters();
            i++;
#endif
            //update reaction rates
            mins = chrono::high_resolution_clock::now();
#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            mine= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_rxn3(mine - mins);
            rxnratetime += elapsed_rxn3.count();
#ifdef CHEMISTRY
            //activate/deactivate compartments
            mins = chrono::high_resolution_clock::now();
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_runspl(mine - mins);
            specialtime += elapsed_runspl.count();

            // update neighbor lists & Binding Managers
            if(tauLastNeighborList >= _neighborListTime) {
                mins = chrono::high_resolution_clock::now();
                rearrangeAllDatabases();
                updateNeighborLists();
                tauLastNeighborList = 0.0;
                mine= chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_runnl2(mine - mins);
                nltime += elapsed_runnl2.count();
            }
            //Special protocols
            mins = chrono::high_resolution_clock::now();
            //special protocols
            executeSpecialProtocols();
            mine= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_runspl2(mine - mins);
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
#endif
        }
    }
    //if run steps were specified, use this
    if(_runSteps != 0) {

#ifdef CHEMISTRY
        while(totalSteps <= _runSteps) {
            //run ccontroller
            if(!_cController->runSteps(_minimizationSteps)) {
                for(auto o: _outputs) o->print(i);
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
                _mController->run();

                // Membrane remeshing
                membraneAdaptiveRemesh();

                updatePositions();
                
                // Update activation of the compartments
                updateActiveCompartments();

                stepsLastMinimization = 0;
            }

            if(stepsLastSnapshot >= _snapshotSteps) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                resetCounters();
                i++;
                stepsLastSnapshot = 0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            resetCounters();
            i++;
#endif

#ifdef DYNAMICRATES
            updateReactionRates();
#endif

#ifdef CHEMISTRY
            //activate/deactivate compartments
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);

            // update neighbor lists
            if(stepsLastNeighborList >= _neighborListSteps) {
                rearrangeAllDatabases();
                updateNeighborLists();
                stepsLastNeighborList = 0;
            }



            //special protocols
            executeSpecialProtocols();
#endif
        }
    }

    //print last snapshots
    for(auto o: _outputs) o->print(i);
    resetCounters();

    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    cout<< "Chemistry time for run=" << chemistrytime <<endl;
    cout << "Minimization time for run=" << minimizationtime <<endl;
    cout<< "Neighbor list + Bmgr time for run="<<nltime<<endl;
    cout<< "Neighbor list time "<<nl2time<<endl;
    cout<< "Bmgr vec time "<<bmgrvectime<<endl;
    cout<< "HYBD time "<<SubSystem::HYBDtime<<endl;
    cout<< "Bmgr time "<<bmgrtime<<endl;
    cout<<"rxnrate time for run="<<rxnratetime<<endl;
    cout<<"Output time for run="<<outputtime<<endl;
    cout<<"Special time for run="<<specialtime<<endl;
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
#ifdef CUDAACCL
    cudaDeviceReset();
#endif
}
