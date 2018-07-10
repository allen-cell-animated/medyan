
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
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
#include "BranchingPoint.h"
#include "Bubble.h"
#include "MTOC.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "MController.h"
#include "Cylinder.h"
#include <unordered_map>
#include 	<tuple>
#include <vector>
#include <algorithm>
#include "Restart.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
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
    _outputs.push_back(new ReactionOut(_outputDirectory + "monomers.traj", _subSystem));
    //Qin add br force out and local diffussing species concentration
    _outputs.push_back(new BRForces(_outputDirectory + "repulsion.traj", _subSystem));
    //_outputs.push_back(new PinForces(_outputDirectory + "pinforce.traj", _subSystem));

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
    
    //Read filament setup, parse filament input file if needed
    FilamentSetup FSetup = p.readFilamentSetup();
//    FilamentData filaments;
    
    cout << "---" << endl;
    cout << "Initializing filaments...";
    
    if(FSetup.inputFile != "") {
        FilamentParser fp(_inputDirectory + FSetup.inputFile);
        filaments = fp.readFilaments();
    }
    fil=get<0>(filaments);
    //add other filaments if specified
    FilamentInitializer* fInit = new RandomFilamentDist();
    
    auto filamentsGen = fInit->createFilaments(_subSystem->getBoundary(),
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

        if(SysParams::Mechanics().transfershareaxis>=0){
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
    //    nvtxRangePushA("neighborlist");
    //Full reset of neighbor lists
    _subSystem->resetNeighborLists();
//    nvtxRangePop();
#ifdef CHEMISTRY
//    nvtxRangePushA("bindingmanager");
    _subSystem->updateBindingManagers();
//    nvtxRangePop();
#endif
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

        mins = chrono::high_resolution_clock::now();
        cout<<"Minimizing energy"<<endl;
//        nvtxRangePushA("mechanics_i");
        _mController->run(false);
//        nvtxRangePop();
        mine= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runm(mine - mins);
        minimizationtime += elapsed_runm.count();
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
#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    oldTau = 0;
#endif

#ifdef MECHANICS
    cout<<"Minimizing energy"<<endl;
    mins = chrono::high_resolution_clock::now();
//    nvtxRangePushA("mechanics_i2");
    _mController->run(false);
//    nvtxRangePop();
    mine= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runm2(mine - mins);
    minimizationtime += elapsed_runm2.count();
//    nvtxRangePushA("output");
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
//    nvtxRangePop();

    cout << "Starting simulation..." << endl;

    int i = 1;

    //if runtime was specified, use this
    if(!areEqual(_runTime, 0.0)) {

#ifdef CHEMISTRY
        //activate/deactivate compartments
        activatedeactivateComp();
        while(tau() <= _runTime) {
            //run ccontroller
//            nvtxRangePushA("chemistry");
            auto var = !_cController->run(_minimizationTime);
//            nvtxRangePop();
            if(var) {
//                nvtxRangePushA("output");
                for(auto o: _outputs) o->print(i);
//                nvtxRangePop();
                break;
            }
            
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
//            std::cout<<endl;
//            std::cout<<"TIME "<<tau()<<endl;
//            std::cout<<endl;
            if(tauLastMinimization >= _minimizationTime) {
                mins = chrono::high_resolution_clock::now();
//                nvtxRangePushA("mechanics");
                _mController->run();
//                nvtxRangePop();
                mine= chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_runm3(mine - mins);
                minimizationtime += elapsed_runm3.count();
//                nvtxRangePushA("update_pos");
                updatePositions();
//                nvtxRangePop();
                tauLastMinimization = 0.0;
            }

            if(tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
//                nvtxRangePushA("output");
                for(auto o: _outputs) o->print(i);
//                nvtxRangePop();
                i++;
                tauLastSnapshot = 0.0;
            }
#elif defined(MECHANICS)
            nvtxRangePushA("output");
            for(auto o: _outputs) o->print(i);
            i++;
            nvtxRangePop();
#endif

#ifdef DYNAMICRATES
//            nvtxRangePushA("rate");
            updateReactionRates();
//            nvtxRangePop();
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(tauLastNeighborList >= _neighborListTime) {
//                nvtxRangePushA("NL");
                updateNeighborLists();
//                nvtxRangePop();
                tauLastNeighborList = 0.0;
            }


            //activate/deactivate compartments
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);

            //special protocols
            executeSpecialProtocols();

            oldTau = tau();
#ifdef CUDAACCL
//            nvtxRangePushA("device reset");
            //reset CUDA context
//            CUDAcommon::handleerror(cudaDeviceSynchronize(), "cudaDeviceSynchronize", "Controller.cu");

//            CUDAcommon::handleerror(cudaDeviceReset(), "cudaDeviceReset", "Controller.cu");

//            nvtxRangePop();
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
            if(!_cController->runSteps(_minimizationSteps)) {
                for(auto o: _outputs) o->print(i);
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
                
                stepsLastMinimization = 0;
            }
            
            if(stepsLastSnapshot >= _snapshotSteps) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                i++;
                stepsLastSnapshot = 0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            i++;
#endif
            
#ifdef DYNAMICRATES
            updateReactionRates();
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
        }
#endif
    }
    
    //print last snapshots
//    nvtxRangePushA("output");
    for(auto o: _outputs) o->print(i);
//    nvtxRangePop();

    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(chk2-chk1);
    cout << "Minimization time for run=" << minimizationtime <<endl;
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout << "Done with simulation!" << endl;
#ifdef CUDAACCL
    cudaDeviceReset();
#endif
}

