
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

#ifndef MEDYAN_ChemCallbacks_h
#define MEDYAN_ChemCallbacks_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Bubble.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Boundary.h"

#include "BindingManager.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"
#include "DissipationTracker.h"

#include <chrono>
#include "CUDAcommon.h"

using namespace mathfunc;

#ifdef RSPECIES_SIGNALING

/// @note - A NOTE TO ALL DEVELOPERS:
///
///         When creating a RSpecies or ReactionBase callback, be sure to not include
///         in the callback structure any variables that may change in the duration
///         of the simulation (i.e. compartments, binding managers, etc). This information
///         needs to be found dynamically to ensure correct behavior - the alternative may
///         cause brutal bugs that are difficult to track in simulation.
///

/// Callback to update the compartment-local binding species based on
/// a change of copy number for an empty site.
struct UpdateBrancherBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateBrancherBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<BranchingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }
                
                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
		CUDAcommon::ctime.tUpdateBrancherBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateBrancherBindingCallback++;
    }
};

struct UpdateLinkerBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateLinkerBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<LinkerBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }

                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateLinkerBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateLinkerBindingCallback++;
    }
};

struct UpdateMotorBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateMotorBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<MotorBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }

                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#ifdef NLSTENCILLIST
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateMotorBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateMotorBindingCallback++;
    }
};

struct UpdateMotorIDCallback{
    
    int _motorType; ///< type of motor to find its binding manager

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateMotorIDCallback(int motorType) : _motorType(motorType) {};
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();
    //DEPRECATED AS OF 9/8/16
//        //find compartment and binding manager
//        Compartment *c = static_cast<Compartment*>(r->getSpecies().getParent());
//        MotorBindingManager* mManager = c->getMotorBindingManager(_motorType);
//        
//        if(delta == +1) {
//            //pull from motorDB of transferred ID's
//            //note that if there is no transfer ID, we are experiencing an unbinding event. The specific
//            //motor will give its ID to the corresponding binding manager.
//            int ID = MotorGhost::_motorGhosts.getTransferID();
//            
//            if(ID != -1) {
//                mManager->addUnboundID(ID);
//                
//                //we can only check this assertion of it is a diffusion event
//                assert(r->getN() == mManager->getAllUnboundIDs().size() &&
//                       "Major bug: number of unbound ID's and copy number does not match");
//            }
//            //else - create an ID. This is an addition at runtime
//            else{
//                mManager->addUnboundID(MotorGhost::_motorGhosts.getId());
//            }
//        }
//        
//        else{ /* -1 */
//            //add to the motorDB of transferred ID's
//            
//            int ID = mManager->getUnboundID();
//            MotorGhost::_motorGhosts.setTransferID(ID);
//            
//            assert(r->getN() == mManager->getAllUnboundIDs().size() &&
//                   "Major bug: number of unbound ID's and copy number does not match");
//        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateMotorIDCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateMotorIDCallback++;
    }
};


#endif

#ifdef REACTION_SIGNALING

/// Callback to extend the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionPlusEndCallback {
    
    Cylinder* _cylinder;
    
    short _plusEnd; ///< Plus end species to mark

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentExtensionPlusEndCallback(Cylinder* cylinder, short plusEnd)
    : _cylinder(cylinder), _plusEnd(plusEnd){};
    
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        //extend the front
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendPlusEnd(_plusEnd);
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentExtensionPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentExtensionPlusEndCallback++;
    }
};

/// Callback to extend the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionMinusEndCallback {
    
    Cylinder* _cylinder;
    
    short _minusEnd; ///< Minus end species to mark

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentExtensionMinusEndCallback(Cylinder* cylinder, short minusEnd)
    : _cylinder(cylinder), _minusEnd(minusEnd){};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        //extend the back
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendMinusEnd(_minusEnd);
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentExtensionMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentExtensionMinusEndCallback++;
    }
};

/// Callback to retract the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionPlusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentRetractionPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f =(Filament*)( _cylinder->getParent());
        f->retractPlusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentRetractionPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentRetractionPlusEndCallback++;
    }
};

/// Callback to retract the minus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionMinusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentRetractionMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->retractMinusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentRetractionMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentRetractionMinusEndCallback++;
    }
};

/// Callback to polymerize the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentPolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizePlusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentPolymerizationPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentPolymerizationPlusEndCallback++;
    }
};

/// Callback to polymerize the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentPolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizeMinusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentPolymerizationMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentPolymerizationMinusEndCallback++;
    }
};

/// Callback to depolymerize the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentDepolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        
        f->depolymerizePlusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDepolymerizationPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDepolymerizationPlusEndCallback++;
    }
};

/// Callback to depolymerize the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentDepolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->depolymerizeMinusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDepolymerizationMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDepolymerizationMinusEndCallback++;
    }
};

/// Callback to unbind a BranchingPoint from a Filament
struct BranchingPointUnbindingCallback {
    
    SubSystem* _ps;
    BranchingPoint* _branchingPoint;

    chrono::high_resolution_clock::time_point mins, mine;

    BranchingPointUnbindingCallback(BranchingPoint* b, SubSystem* ps)
    : _ps(ps), _branchingPoint(b) {}
    
    void operator() (ReactionBase *r) {
	    mins = chrono::high_resolution_clock::now();

//        //@{
//        std::cout<<"Brancher unbinding "<<_branchingPoint->getFirstCylinder()->getId()<<" "
//                ""<<_branchingPoint->getPosition()
//                 <<" "<<_branchingPoint->getSecondCylinder()->getId()<<endl;
//        //@}
        //remove the branching point
        _ps->removeTrackable<BranchingPoint>(_branchingPoint);
        delete _branchingPoint;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tBranchingPointUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cBranchingPointUnbindingCallback++;
    }
    
    
};


/// Callback to create a BranchingPoint on a Filament
struct BranchingCallback {
    
    SubSystem* _ps;        ///< ptr to subsystem
    
    BranchingManager* _bManager; ///< Branching manager for this compartment
    
    short _plusEnd;        ///< Plus end marker of new cylinder
    
    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction

    chrono::high_resolution_clock::time_point mins, mine;

    BranchingCallback(BranchingManager* bManager, short plusEnd,
                      float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _bManager(bManager),
    _plusEnd(plusEnd), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
	    mins = chrono::high_resolution_clock::now();
        BranchingPoint* b;
        float frate;
        short branchType = _bManager->getBoundInt();
        
        //choose a random binding site from manager
        tuple<CCylinder*, short> site;
#ifdef NLORIGINAL
        site = _bManager->chooseBindingSite();
#endif
#ifdef NLSTENCILLIST
        site = _bManager->chooseBindingSitesstencil();
#endif
        
        //get info from site
        Cylinder* c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();
        
        floatingpoint pos = floatingpoint(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];
        if(SysParams::RUNSTATE==true){
        //Get a position and direction of a new filament
        auto x1 = c1->getFirstBead()->vcoordinate();
        auto x2 = c1->getSecondBead()->vcoordinate();
        
        //get original direction of cylinder
        auto p= midPointCoordinate(x1, x2, pos);
        vector<floatingpoint> n = twoPointDirection(x1, x2);
        
        //get branch projection
#ifdef MECHANICS
        //use mechanical parameters
        floatingpoint l, t;
        if(SysParams::Mechanics().BrStretchingL.size() != 0) {
            l = SysParams::Mechanics().BrStretchingL[branchType];
            t = SysParams::Mechanics().BrBendingTheta[branchType];
        }
        else {
            cout << "Branching initialization cannot occur unless mechanical parameters are specified."
            << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
            << endl;
            l = 10.0;
            t = 1.22;
        }
#else
        cout << "Branching initialization cannot occur unless mechanics is enabled. Using"
        << " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
        << endl;
        floatingpoint l = 10.0;
        floatingpoint t = 1.22;
#endif
        floatingpoint s = SysParams::Geometry().monomerSize[filType];
        
        auto branchPosDir = branchProjection(n, p, l, s, t);
        auto bd = get<0>(branchPosDir); auto bp = get<1>(branchPosDir);
        
        //create a new filament
        Filament* f = _ps->addTrackable<Filament>(_ps, filType, bp, bd, true, true);
        
        //mark first cylinder
        Cylinder* c = f->getCylinderVector().front();
        c->getCCylinder()->getCMonomer(0)->speciesPlusEnd(_plusEnd)->up();
        
        //create new branch
            
        b= _ps->addTrackable<BranchingPoint>(c1, c, branchType, pos);
            
        frate=_offRate;
        }
        else
        {
            CCylinder* c = nullptr; auto check = false;
            vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> BrT=_bManager->getbtuple();
            for(auto T:BrT){
                CCylinder* cx=get<0>(get<0>(T));
                floatingpoint p = floatingpoint(get<1>(get<0>(T)))/ floatingpoint(SysParams::Geometry().cylinderNumMon[filType]);
                if(cx->getCylinder()->getId()==c1->getId() && p==pos){
                    c=get<0>(get<1>(T));
                    check = true;
                    break;
                }}
            if(check){

                
            b= _ps->addTrackable<BranchingPoint>(c1, c->getCylinder(), branchType, pos);
                
            frate=0.0;
            }
            else {
                b = nullptr;
                cout<<"Brancher Error. Cannot find binding Site in the list. Cannot complete restart. Exiting." <<endl;
                exit(EXIT_FAILURE);
            }
        }
        
        //create off reaction
        auto cBrancher = b->getCBranchingPoint();
        
        cBrancher->setRates(_onRate, frate);
        cBrancher->createOffReaction(r, _ps);
        cBrancher->getOffReaction()->setBareRate(SysParams::BUBBareRate[branchType]);
#ifdef DYNAMICRATES

        b -> updateReactionRates();
#endif
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tBranchingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cBranchingCallback++;
    }
};



/// Callback to unbind a Linker from a Filament
struct LinkerUnbindingCallback {
    
    SubSystem* _ps;
    Linker* _linker;

    chrono::high_resolution_clock::time_point mins, mine;

    LinkerUnbindingCallback(Linker* l, SubSystem* ps) : _ps(ps), _linker(l) {}
    
    void operator() (ReactionBase *r) {
	    CUDAcommon::tmin.linkerunbindingcalls++;
	    mins = chrono::high_resolution_clock::now();

//        //@{
//        std::cout<<"Linker unbinding "<<_linker->getFirstCylinder()->getId()<<" "<<_linker->getFirstPosition()
//                 <<" "<<_linker->getSecondCylinder()->getId()<<" "<<_linker->getSecondPosition()<<endl;
//        //@}
        //remove the linker
        _ps->removeTrackable<Linker>(_linker);
        delete _linker;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tLinkerUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cLinkerUnbindingCallback++;
    }
};

/// Callback to bind a Linker to Filament
struct LinkerBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    LinkerBindingManager* _lManager; ///< Linker binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    DissipationTracker* _dt;

    chrono::high_resolution_clock::time_point mins, mine;

    LinkerBindingCallback(LinkerBindingManager* lManager,
                          float onRate, float offRate, SubSystem* ps, DissipationTracker* dt)
    
    : _ps(ps), _lManager(lManager), _onRate(onRate), _offRate(offRate), _dt(dt) {}
    
    void operator() (ReactionBase *r) {
	    CUDAcommon::tmin.linkerbindingcalls++;
	    mins = chrono::high_resolution_clock::now();

        //get a random binding
        short linkerType = _lManager->getBoundInt();
        float f;
        
        //choose a random binding site from manager
        vector<tuple<CCylinder*, short>> site;
#ifdef NLORIGINAL
        site = _lManager->chooseBindingSites();
#endif
#ifdef NLSTENCILLIST
        site = _lManager->chooseBindingSitesstencil();
#endif
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        short filType = c1->getType();
        
        // Create a linker
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        
        floatingpoint pos1 = floatingpoint(get<1>(site[0])) / cylinderSize;
        floatingpoint pos2 = floatingpoint(get<1>(site[1])) / cylinderSize;
        
        Linker* l = _ps->addTrackable<Linker>(c1, c2, linkerType, pos1, pos2);
        
        //create off reaction
        auto cLinker = l->getCLinker();
        //aravind June 24, 2016.
        if(SysParams::RUNSTATE==false)
            f=0.0;
        else
            f=_offRate;
        //@
        cLinker->setRates(_onRate, f);
        cLinker->createOffReaction(r, _ps);
        
        if(SysParams::Chemistry().eventTracking){
            _dt->recordLinkerBinding(l);
        }

#ifdef DYNAMICRATES
        //reset the associated reactions
        //aravind june 24, 2016
        cLinker->getOffReaction()->setBareRate(SysParams::LUBBareRate[linkerType]);
        //@
        l->updateReactionRates();
#endif
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tLinkerBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cLinkerBindingCallback++;
    }
};

/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    MotorGhost* _motor;

    chrono::high_resolution_clock::time_point mins, mine;

    MotorUnbindingCallback(MotorGhost* m, SubSystem* ps) :
    
    _ps(ps), _motor(m) {}
    
    void operator() (ReactionBase *r) {
    	CUDAcommon::tmin.motorunbindingcalls++;
	    mins = chrono::high_resolution_clock::now();

        //find the motor binding manager associated with this species
//        Species* sd = &(r->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
//        Compartment* c = static_cast<Compartment*>(sd->getParent());
//        auto mManager = c->getMotorBindingManager(_motor->getType());
        
        //DEPRECATED AS OF 9/8/16
//        
//        mManager->removeUnboundID(MotorGhost::_motorGhosts.deleteID());
//        
//        //re-add unbound ID to motor binding manager
//        mManager->addUnboundID(_motor->getId());

//        //@{
//        std::cout<<"Motor unbinding "<<_motor->getFirstCylinder()->getId()<<" "<<_motor->getFirstPosition()
//                  <<" "<<_motor->getSecondCylinder()->getId()<<" "<<_motor->getSecondPosition()<<endl;
//        //@}

        //remove the motor
        _ps->removeTrackable<MotorGhost>(_motor);
        delete _motor;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorUnbindingCallback++;
    }
};


/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    MotorBindingManager* _mManager;///< Motor binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    chrono::high_resolution_clock::time_point mins, mine;

    MotorBindingCallback(MotorBindingManager* mManager,
                         float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _mManager(mManager), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
	    CUDAcommon::tmin.motorbindingcalls++;
	    mins = chrono::high_resolution_clock::now();

        //get a random binding
        short motorType = _mManager->getBoundInt();
        float f;
        //choose a random binding site from manager
        vector<tuple<CCylinder*, short>> site;
#ifdef NLORIGINAL
        site = _mManager->chooseBindingSites();
#endif
#ifdef NLSTENCILLIST
        site = _mManager->chooseBindingSitesstencil();
#endif
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        short filType = c1->getType();
        
        // Create a motor
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        
        floatingpoint pos1 = floatingpoint(get<1>(site[0])) / cylinderSize;
        floatingpoint pos2 = floatingpoint(get<1>(site[1])) / cylinderSize;
        
        MotorGhost* m = _ps->addTrackable<MotorGhost>(c1, c2, motorType, pos1, pos2, _onRate, _offRate);
        
        //attach an ID to the motor based on the transfer ID
        //DEPRECATED AS OF 9/22/16
//        m->setID(MotorGhost::_motorGhosts.getTransferID());
        
        //create off reaction
        auto cMotorGhost = m->getCMotorGhost();
        //aravind June 24, 2016.
        if(SysParams::RUNSTATE==false){
        f=0.0;
        }
        else
            f=_offRate;
        //@
        cMotorGhost->setRates(_onRate, f);
        cMotorGhost->createOffReaction(r, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated walking reactions
        m->updateReactionRates();
        //aravind June 24,2016.
        cMotorGhost->getOffReaction()->setBareRate(SysParams::MUBBareRate[motorType]);
        //@

#endif
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorBindingCallback++;

    }
};

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingCallback {
    
    Cylinder* _c;        ///< Cylinder this callback is attached to
    
    short _oldPosition;  ///< Old position of motor head
    short _newPosition;  ///< New position of motor head
    
    short _motorType;    ///< Type of motor
    short _boundType;    ///< Type of bound this motor took place of
    
    SubSystem* _ps;      ///< Ptr to subsystem

    chrono::high_resolution_clock::time_point mins, mine;
    DissipationTracker* _dt;

    MotorWalkingCallback(Cylinder* c,
                         short oldPosition, short newPosition,
                         short motorType, short boundType, SubSystem* ps, DissipationTracker* dt)
    
    :_c(c), _oldPosition(oldPosition), _newPosition(newPosition),
    _motorType(motorType), _boundType(boundType), _ps(ps), _dt(dt) {}
    
    void operator() (ReactionBase* r) {
	    CUDAcommon::tmin.motorwalkingcalls++;
	    mins = chrono::high_resolution_clock::now();
//        cout<<"Motor walking begins"<<endl;
        //get species
        CCylinder* cc = _c->getCCylinder();
        CMonomer* monomer = cc->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        
        short filType = _c->getType();
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();


        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];

        floatingpoint oldpos = floatingpoint(_oldPosition) / cylinderSize;
        floatingpoint newpos = floatingpoint(_newPosition) / cylinderSize;
#ifdef CROSSCHECK
        auto sb = SysParams::Chemistry().motorBoundIndex[0];
	    CMonomer* monomernew = cc->getCMonomer(_newPosition);

	    cout<<"mw before speciesMotor E-0/B-1 "<<sm1->getN()<<" "
	    <<monomernew->speciesMotor(_motorType)->getN()<<endl;
	    cout<<"mw before speciesBound E-1/B-0 "<<monomer->speciesBound(sb)->getN()
	    <<" "<<monomernew->speciesBound(sb)->getN()<<endl;
#endif
	    m->moveMotorHead(_c, oldpos, newpos, _boundType, _ps);
#ifdef CROSSCHECK
	    cout<<"mw after speciesMotor E-0/B-1 "<<sm1->getN()<<" "<<
	    monomernew->speciesMotor(_motorType)->getN()<<endl;
	    cout<<"mw after speciesBound E-1/B-0 "<<monomer->speciesBound(sb)->getN()
	        <<" "<<monomernew->speciesBound(sb)->getN()<<endl;
#endif
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        m->updateReactionRates();

#endif
        
        if(SysParams::Chemistry().eventTracking){
            _dt->recordWalk(m);
        }
        
        
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorWalkingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorWalkingCallback++;
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderCallback {
    
    Cylinder* _oldC;        ///< Old cylinder the motor is attached to
    Cylinder* _newC;        ///< New cylinder motor will be attached to
    
    short _oldPosition;     ///< Old position of motor head
    short _newPosition;     ///< New position of motor head
    
    short _motorType;       ///< Type of motor
    short _boundType;       ///< Type of bound this motor is taking place of
    
    SubSystem* _ps;         ///< Ptr to subsystem

    chrono::high_resolution_clock::time_point mins, mine;

    MotorMovingCylinderCallback(Cylinder* oldC, Cylinder* newC,
                                short oldPosition, short newPosition,
                                short motorType, short boundType, SubSystem* ps)
    
    :_oldC(oldC), _newC(newC), _oldPosition(oldPosition), _newPosition(newPosition),
    _motorType(motorType), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
	    CUDAcommon::tmin.motorwalkingcalls++;
	    mins = chrono::high_resolution_clock::now();
//        cout<<"Motor moving cylinder begins"<<endl;
        //get species
        CCylinder* oldCC = _oldC->getCCylinder();
        CMonomer* monomer = oldCC->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        short filType = _oldC->getType();

        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
/*        cout<<"filament Type "<<filType<<endl;
        cout<<"cylinder size "<<cylinderSize<<endl;
        cout<<"Cylinder oldpos "<<_oldC->getID()<<" "<<_oldPosition<<endl;
        cout<<"Cylinder newpos "<<_newC->getID()<<" "<<_newPosition<<endl;
        cout<<"-----------"<<endl;*/
        floatingpoint oldpos = floatingpoint(_oldPosition) / cylinderSize;
        floatingpoint newpos = floatingpoint(_newPosition) / cylinderSize;
        
        m->moveMotorHead(_oldC, _newC, oldpos, newpos, _boundType, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        m->updateReactionRates();
#endif
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorMovingCylinderCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorMovingCylinderCallback++;
    }
};

/// Struct to create a new filament based on a given reaction
struct FilamentCreationCallback {
    
    //@{
    /// Integer specifying the chemical properties of first cylinder
    short _plusEnd;
    short _minusEnd;
    short _filament;
    //@}
    
    ///Filament type to create
    short _filType;
    
    Compartment* _compartment; ///< compartment to put this filament in
    SubSystem* _ps; ///< Ptr to the subsystem

    chrono::high_resolution_clock::time_point mins, mine;

    FilamentCreationCallback(short plusEnd, short minusEnd, short filament,
                             short filType, SubSystem* ps, Compartment* c = nullptr)
    
    : _plusEnd(plusEnd), _minusEnd(minusEnd), _filament(filament),
    _filType(filType), _compartment(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
	    mins = chrono::high_resolution_clock::now();

        Compartment* c;
        
        //no compartment was set, pick a random one
        if(_compartment == nullptr)
            c = GController::getRandomCompartment();
        else
            c = _compartment;
        
        //set up a random initial position and direction
        vector<floatingpoint> position;
        vector<floatingpoint> direction;

        if(_ps->getBoundary()->getShape() == BoundaryShape::Cylinder) {
            while(true) {
                position = GController::getRandomCenterCoordinates(c);
                
                //getting random numbers between -1 and 1
                
                direction = {Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1), 0};
                normalize(direction);
                
                auto npp = nextPointProjection(position,
                                               SysParams::Geometry().cylinderSize[_filType], direction);
                
                bool inbubble = false;
                //check if within bubble
                for(auto bb : Bubble::getBubbles()) {
                    auto radius = bb->getRadius();

                    if((twoPointDistancesquared(bb->getBead()->vcoordinate(), position) < (radius * radius)) ||
                       (twoPointDistancesquared(bb->getBead()->vcoordinate(), npp) < (radius * radius))){
                        inbubble = true;
                        break;
                    }
                }

                if(inbubble) break;

                //check if within boundary && if within REPULSIONEXPIN region (set as 125nm)
                if(_ps->getBoundary()->within(position) &&
                   _ps->getBoundary()->within(npp)){
                    if(_ps->getBoundary()->distance(position) > 125 &&
                       _ps->getBoundary()->distance(npp) > 125)
                    	break;
                }

            }
            
            //create filament, set up ends and filament species
            Filament* f = _ps->addTrackable<Filament>(_ps, _filType, position, direction, true, false);
            
            //initialize the nucleation
            f->nucleate(_plusEnd, _filament, _minusEnd);
        }
        else {
            while(true) {
                position = GController::getRandomCoordinates(c);
                
                //getting random numbers between -1 and 1
                
                direction = {Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1)};
                normalize(direction);
                
                auto npp = nextPointProjection(position,
                                               SysParams::Geometry().cylinderSize[_filType], direction);
                
                bool inbubble = false;
                //check if within bubble
                for(auto bb : Bubble::getBubbles()) {
                    auto radius = bb->getRadius();

                    if((twoPointDistancesquared(bb->getBead()->vcoordinate(), position) < (radius * radius)) ||
                       (twoPointDistancesquared(bb->getBead()->vcoordinate(), npp) < (radius * radius))){
                        inbubble = true;
                        break;
                    }
                }

                if(inbubble) break;

                //check if within boundary
                if(_ps->getBoundary()->within(position) &&
                   _ps->getBoundary()->within(npp))
                    break;
            }
            
            //create filament, set up ends and filament species
            Filament* f = _ps->addTrackable<Filament>(_ps, _filType, position, direction, true, false);
            
            //initialize the nucleation
            f->nucleate(_plusEnd, _filament, _minusEnd);
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentCreationCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentCreationCallback++;
        }

    
};

///Struct to sever a filament based on a reaction
struct FilamentSeveringCallback {
    
    Cylinder* _c1;  ///< Filament severing point

    chrono::high_resolution_clock::time_point mins, mine;

    FilamentSeveringCallback(Cylinder* c1) : _c1(c1) {}
    
    void operator() (ReactionBase* r) {
	    mins = chrono::high_resolution_clock::now();

        //reactants should be re-marked
        for(int i = 0; i < SEVERINGREACTANTS + 1; i++)
            r->rspecies()[i]->up();
        
        //sever the filament at given position
        Filament* f = (Filament*)(_c1->getParent());
        
        //create a new filament by severing
        Filament* newF = f->sever(_c1->getPosition());
        
        //if severing did not occur, return
        if(newF == nullptr) return;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentSeveringCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentSeveringCallback++;
    }
};

/// Struct to destroy a filament based on a reaction
struct FilamentDestructionCallback {
    
    Cylinder* _c; ///< Cylinder to destroy
    
    SubSystem* _ps; ///< SubSystem ptr

    chrono::high_resolution_clock::time_point mins, mine;

    FilamentDestructionCallback(Cylinder* c, SubSystem* ps) : _c(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
	    mins = chrono::high_resolution_clock::now();

        Filament* f = (Filament*)(_c->getParent());
        
        _ps->removeTrackable<Filament>(f);
        delete f;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDestructionCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDestructionCallback++;
    }
};

#endif

#endif
