
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

#include "core/controller/Controller.h"

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
#include "Membrane.hpp"
#include "MembraneHierarchy.h"
#include "MembraneRegion.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "core/controller/MController.h"
#include "Cylinder.h"
#include <unordered_map>
#include <tuple>
#include <vector>
#include <algorithm>
#include "Rand.h"
#include "Restart.h"
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
    
    //general check of macros
#if defined(DYNAMICRATES) && (!defined(CHEMISTRY) || !defined(MECHANICS))
    cout << "If dynamic rates is turned on, chemistry and mechanics must be "
         << "defined. Please set these compilation macros and try again. Exiting."
         << endl;
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

    //Always read geometry, check consistency
    p.readGeoParams();
    if(!SysParams::checkGeoParameters()) exit(EXIT_FAILURE);
    
    //CALLING ALL CONTROLLERS TO INITIALIZE
    //Initialize geometry controller
    cout << "---" << endl;
    cout << "Initializing geometry...";
    _gController->initializeGrid();
    cout << "Done." << endl;
    
    //Initialize boundary
    cout << "---" << endl;
    cout << "Initializing boundary...";
    
    auto BTypes = p.readBoundaryType();
    p.readBoundParams();
    
    //initialize
    _gController->initializeBoundary(BTypes);
    cout << "Done." <<endl;
    
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
        cout << "Need to specify a chemical input file. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    _cController->initialize(CAlgorithm.algorithm, ChemData);
    cout << "Done." << endl;
    
    //Set up chemistry output if any
    string chemsnapname = _outputDirectory + "chemistry.traj";
    _outputs.push_back(new Chemistry(chemsnapname, _subSystem, ChemData,
                                     _subSystem->getCompartmentGrid()));

    string concenname = _outputDirectory + "concentration.traj";
    _outputs.push_back(new Concentrations(concenname, _subSystem, ChemData));

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

    /**************************************************************************
    Now starting to add the membrane into the network.
    **************************************************************************/
    MembraneSetup MemSetup = p.readMembraneSetup();
    
    cout << "---" << endl;
    cout << "Initializing membranes...";
    
    if(MemSetup.inputFile != "") {
        MembraneParser memp(_inputDirectory + MemSetup.inputFile);
        membraneData = memp.readMembranes();
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

        Membrane* newMembrane = _subSystem->addTrackable<Membrane>(_subSystem, type, it);

        // Update membrane geometry
        newMembrane->updateGeometry(true);

        // Add to the global membrane hierarchy
        MembraneHierarchy::addMembrane(newMembrane);
    }
    cout << "Done. " << membraneData.size() << " membranes created." << endl;

    // Create a region inside the membrane
    auto regionInMembrane = (
        membraneData.size()?
        MembraneRegion::makeByChildren(&MembraneHierarchy::getRoot()):
        make_unique<MembraneRegion>(_subSystem->getBoundary())
    );

    // Deactivate all the compartments outside membrane, and mark boundaries as interesting
    for(auto c : _subSystem->getCompartmentGrid()->getCompartments()) {
        if(!c->getTriangles().empty()) {
            // Contains triangles, so this compartment is at the boundary.
            c->boundaryInteresting = true;

            // Update partial activate status
            c->getSlicedVolumeArea();
            _cController->updateActivation(c);

        } else if( ! regionInMembrane->contains(vector2Array<double, 3>(c->coordinates()))) {
            // Compartment is outside the membrane
            _cController->deactivate(c);
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
                        auto tc = ac[Rand::randInteger(1, nac-1)];
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
        int numSegment = d / SysParams::Geometry().cylinderSize[type];
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

        if(SysParams::Mechanics().transfershareaxis!=3){
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

        }
}
void Controller::ControlfrontbackEndComp(){
    Compartment* maxcomp=NULL;
    Compartment* mincomp=NULL;
//    Bead* minbead=NULL;

    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
        auto cyls=C->getCylinders();
        if(cyls.size()>0){
            if(maxcomp==NULL)
                maxcomp=C;
            else{
                auto mcoord=maxcomp->coordinates();
                auto ccord=C->coordinates();
                if(mcoord[SysParams::Mechanics().transfershareaxis]<ccord[SysParams::Mechanics().transfershareaxis])
                    maxcomp=C;
            }
            if(mincomp==NULL)
                mincomp=C;
            else{
                auto mcoord=mincomp->coordinates();
                auto ccord=C->coordinates();
                if(mcoord[SysParams::Mechanics().transfershareaxis]>ccord[SysParams::Mechanics().transfershareaxis])
                    mincomp=C;
            }
        }
    }

    // front end
    auto cmaxcomp=maxcomp->coordinates();
    for(auto C:maxcomp->getNeighbours()){
        auto cC=C->coordinates();
        if(cmaxcomp[SysParams::Mechanics().transfershareaxis]<cC[SysParams::Mechanics().transfershareaxis])
            maxcomp=C;
    }
    cmaxcomp=maxcomp->coordinates();
    for(auto C:maxcomp->getNeighbours()){
        auto cC=C->coordinates();
        if(cmaxcomp[SysParams::Mechanics().transfershareaxis]<cC[SysParams::Mechanics().transfershareaxis])
            maxcomp=C;
    }
    cmaxcomp=maxcomp->coordinates();
    assert((maxcomp!=NULL) && "Non existent maxcomp. Exiting.");
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
        auto cC=C->coordinates();

        if(cC[SysParams::Mechanics().transfershareaxis]>cmaxcomp[SysParams::Mechanics().transfershareaxis]){
            if(C->isActivated())
                fCompmap.insert(pair<int,Compartment*>(cC[SysParams::Mechanics().transfershareaxis],C));
        }
        else{
            if(!(C->isActivated()))
                activatecompartments.push_back(C);
        }
    }
    //back end
    auto cmincomp=mincomp->coordinates();
    for(auto C:mincomp->getNeighbours()){
      auto cC=C->coordinates();
        if(cmincomp[SysParams::Mechanics().transfershareaxis]>cC[SysParams::Mechanics().transfershareaxis])
            mincomp=C;
    }
    cmincomp=mincomp->coordinates();
    for(auto C:mincomp->getNeighbours()){
        auto cC=C->coordinates();
        if(cmincomp[SysParams::Mechanics().transfershareaxis]>cC[SysParams::Mechanics().transfershareaxis])
            mincomp=C;
    }
    cmincomp=mincomp->coordinates();
    assert(mincomp!=NULL && "Non existent mincomp. Exiting.");
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
        auto cC=C->coordinates();
        if(cC[SysParams::Mechanics().transfershareaxis]<cmincomp[SysParams::Mechanics().transfershareaxis]){
            auto it = std::find(activatecompartments.begin(), activatecompartments.end(), C);
            if(it!=activatecompartments.end())
                activatecompartments.erase(it);
            if(C->isActivated()){
                bCompmap.insert(pair<int,Compartment*>(cC[SysParams::Mechanics().transfershareaxis],C));
            }
        }
    }
    std::cout<<"Maxcomp "<<maxcomp->coordinates()[SysParams::Mechanics().transfershareaxis]<<" ";
    std::cout<<"Mincomp "<<mincomp->coordinates()[SysParams::Mechanics().transfershareaxis]<<endl;

}

void Controller::moveBoundary(double deltaTau) {
    
    //calculate distance to move
    double dist = SysParams::Boundaries().moveSpeed * deltaTau;
    if(abs(dist)>0){
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
                c->getSlicedVolumeArea();
                _cController->updateActivation(c);

                // No matter whether the compartment is interesting before, mark it as interesting
                c->boundaryInteresting = true;
            } else if(c->boundaryInteresting) { // Interesting last round but now empty
                bool inMembrane = (
                    (!theMembrane->isClosed()) ||
                    (theMembrane->signedDistance(vector2Array<double, 3>(c->coordinates()), false) < 0.0)
                );
                if(inMembrane) {
                    // Fully activate the compartment
                    c->setPartialVolume(GController::getCompartmentVolume());
					auto& fullArea = GController::getCompartmentArea();
                    c->setPartialArea({{
                        fullArea[0], fullArea[0],
						fullArea[1], fullArea[1],
						fullArea[2], fullArea[2]
                    }});
                    _cController->updateActivation(c);
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
    
    //Full reset of neighbor lists
    _subSystem->resetNeighborLists();
    
#ifdef CHEMISTRY
    _subSystem->updateBindingManagers();
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
            
            cout << _subSystem->getBoundary()->distance(b->coordinate) << endl;
            cout << SysParams::Mechanics().pinDistance << endl;
            
            
            //if within dist to boundary, add
            if(_subSystem->getBoundary()->distance(b->coordinate) < SysParams::Mechanics().pinDistance) {
                
                b->pinnedPosition = b->coordinate;
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

            //cout << _subSystem->getBoundary()->lowerdistance(b->coordinate) << endl;
            //cout << SysParams::Mechanics().pinDistance << endl;

            auto index = Rand::randDouble(0,1);
            //cout << index <<endl;
            //if within dist to boundary and index > 0.5, add
            if(_subSystem->getBoundary()->lowerdistance(b->coordinate) < SysParams::Mechanics().pinDistance
                && index < SysParams::Mechanics().pinFraction && b->isPinned() == false) {
                //cout << index << endl;
                b->pinnedPosition = b->coordinate;
                b->addAsPinned();
            }
        }
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
    chrono::high_resolution_clock::time_point chk1, chk2;
    chk1 = chrono::high_resolution_clock::now();
//RESTART PHASE BEGINS
    if(SysParams::RUNSTATE==false){
        cout<<"RESTART PHASE BEINGS."<<endl;
        Restart* _restart = new Restart(_subSystem, filaments,_chemData);
//Step 1. Turn off diffusion, passivate filament reactions and empty binding managers.
//        _restart->settorestartphase();
        cout<<"Turned off Diffusion, filament reactions."<<endl;
//Step 2. Add bound species to their respective binding managers. Turn off unbinding, update propensities.
        _restart->addtoHeapbranchers();
        _restart->addtoHeaplinkermotor();
        cout<<"Bound species added to reaction heap."<<endl;
//Step 2A. Turn off diffusion, passivate filament reactions and empty binding managers.
                _restart->settorestartphase();
//Step 3. ############ RUN LINKER/MOTOR REACTIONS TO BIND BRANCHERS, LINKERS, MOTORS AT RESPECTIVE POSITIONS.#######
        std::cout<<"Reactions to be fired "<<_restart->getnumchemsteps()<<endl;
        _cController->runSteps(_restart->getnumchemsteps());
        cout<<"Reactions fired! Displaying heap"<<endl;
//Step 4. Display the number of reactions yet to be fired. Should be zero.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto &Mgr:C->getFilamentBindingManagers()){
                if(Mgr->numBindingSites()==0)
                    cout<< Mgr->numBindingSites()<<' ';
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
        
        cout<<"Minimizing energy"<<endl;
        _mController->run(false);
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
#ifdef MECHANICS
     cout<<"Minimizing energy"<<endl;
    _mController->run(false);
    
    //reupdate positions and neighbor lists
    updatePositions();
    updateNeighborLists();
    
#ifdef DYNAMICRATES
    updateReactionRates();
#endif
    
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
        activatedeactivateComp();
        while(tau() <= _runTime) {
            //run ccontroller
            if(!_cController->run(_minimizationTime)) {
                for(auto o: _outputs) o->print(i);
                resetCounters();
                break;
            }
            
            //add the last step
            tauLastSnapshot += tau() - oldTau;
            tauLastMinimization += tau() - oldTau;
            tauLastNeighborList += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            if(tauLastMinimization >= _minimizationTime) {
                _mController->run();
                updatePositions();

                // Update activation of the compartments
                updateActiveCompartments();

#ifdef DYNAMICRATES
                updateReactionRates();
#endif

                tauLastMinimization = 0.0;
            }
            
            if(tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                resetCounters();
                i++;
                tauLastSnapshot = 0.0;
            }

#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            resetCounters();
            i++;
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(tauLastNeighborList >= _neighborListTime) {
                updateNeighborLists();
                tauLastNeighborList = 0.0;
            }
            
            //activate/deactivate compartments
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);

            //special protocols
            executeSpecialProtocols();
            
            oldTau = tau();
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
                updatePositions();
                
                // Update activation of the compartments
                updateActiveCompartments();

#ifdef DYNAMICRATES
                updateReactionRates();
#endif
                
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
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(stepsLastNeighborList >= _neighborListSteps) {
                updateNeighborLists();
                stepsLastNeighborList = 0;
            }
            
            //move the boundary
            moveBoundary(tau() - oldTau);
            //activate/deactivate compartments
            activatedeactivateComp();
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
    
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
}

