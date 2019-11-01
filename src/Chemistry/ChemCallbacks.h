
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

#ifndef MEDYAN_ChemCallbacks_h
#define MEDYAN_ChemCallbacks_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "CaMKIIingPoint.h"
#include "Boundary.h"

#include "BindingManager.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

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
    
    //Constructor, sets members
    UpdateBrancherBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
        
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<BranchingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);
                
                else /* -1 */manager->removePossibleBindings(cc, _bindingSite);
            }
        }
    }
};

struct UpdateCaMKIIerBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update
    
    //Constructor, sets members
    UpdateCaMKIIerBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
        //cout<<"========== CaMKIIBindingManager CallBack" <<endl; //Carlos verbose prints
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();

        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<CaMKIIBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);
                
                else /* -1 */manager->removePossibleBindings(cc, _bindingSite);
            }
        }
    }
};

struct UpdateCaMKIIerBundlingCallback {

    Cylinder* _cylinder; ///< cylinder to update

    short _bindingSite;  ///< binding site to update

    //Constructor, sets members
    UpdateCaMKIIerBundlingCallback(Cylinder* cylinder, short bindingSite)

    : _cylinder(cylinder), _bindingSite(bindingSite) {}

    //callback
    void operator() (RSpecies *r, int delta) {
    	if(delta != 1 && delta != -1)
        	cerr <<"========== CaMKIIBundlingManager CallBack ---- delta = " << delta <<endl; //Carlos verbose prints

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();

        for(auto &manager : c->getFilamentBindingManagers()) {

            if(dynamic_cast<CaMKIIBundlingManager*>(manager.get())) {

                CCylinder* cc = _cylinder->getCCylinder();

                //update binding sites
                if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);

                else /* -1 */ {
                    manager->removePossibleBindings(cc, _bindingSite);

                    cout << "---------------============----------" << endl;
                    for(auto ccc : Cylinder::getCylinders()) {
                        if(dynamic_cast<CaMKIICylinder*>(ccc)) {
                            cout << "========MILLAD: " << ccc->getCCylinder() << endl;
                        }
                    }
                    cout << "---------------============----------" << endl;
                }
            }
        }
    }
};

struct UpdateCaMKIIerDummyCylinderCallback {

    Cylinder* _cylinder; ///< cylinder to update

    short _bindingSite;  ///< binding site to update

    //Constructor, sets members
    UpdateCaMKIIerDummyCylinderCallback(Cylinder* cylinder, short bindingSite)

            : _cylinder(cylinder), _bindingSite(bindingSite) {}

    //callback
    void operator() (RSpecies *r, int delta) {
        if(!dynamic_cast<CaMKIICylinder*>(_cylinder))
            return;

        const int n = r->getN();

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();

        for(auto &manager : c->getFilamentBindingManagers()) {

            if(dynamic_cast<CaMKIIBundlingManager*>(manager.get())) {

                CCylinder* cc = _cylinder->getCCylinder();

                if(delta == -1 && n == 0) {
                    manager->removePossibleBindings(cc, _bindingSite);
                    cerr << "====== MILLAD: Removing CaMKIICylinder from PossibleCylinders " << _cylinder << "  " << _cylinder->getCCylinder()
                         << endl;
                }
            }
        }
    }
};

struct UpdateLinkerBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update
    
    //Constructor, sets members
    UpdateLinkerBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
        
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<LinkerBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);
                
                else /* -1 */manager->removePossibleBindings(cc, _bindingSite);
            }
        }
    }
};

struct UpdateMotorBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update
    
    //Constructor, sets members
    UpdateMotorBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
        
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<MotorBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);
                
                else /* -1 */ manager->removePossibleBindings(cc, _bindingSite);
            }
        }
    }
};

struct UpdateMotorIDCallback{
    
    int _motorType; ///< type of motor to find its binding manager
    
    //Constructor, sets members
    UpdateMotorIDCallback(int motorType) : _motorType(motorType) {};
    
    //callback
    void operator() (RSpecies *r, int delta) {
        
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
//                mManager->addUnboundID(MotorGhost::_motorGhosts.getID());
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
    }
};


#endif

#ifdef REACTION_SIGNALING

/// Callback to extend the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionPlusEndCallback {
    
    Cylinder* _cylinder;
    
    short _plusEnd; ///< Plus end species to mark
    
    //Constructor, sets members
    FilamentExtensionPlusEndCallback(Cylinder* cylinder, short plusEnd)
    : _cylinder(cylinder), _plusEnd(plusEnd){};
    
    //Callback
    void operator() (ReactionBase *r){
        //extend the front
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendPlusEnd(_plusEnd);
    }
};

/// Callback to extend the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionMinusEndCallback {
    
    Cylinder* _cylinder;
    
    short _minusEnd; ///< Minus end species to mark
    
    //Constructor, sets members
    FilamentExtensionMinusEndCallback(Cylinder* cylinder, short minusEnd)
    : _cylinder(cylinder), _minusEnd(minusEnd){};
    //Callback
    void operator() (ReactionBase *r){
        //extend the back
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendMinusEnd(_minusEnd);
    }
};

/// Callback to retract the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionPlusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentRetractionPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f =(Filament*)( _cylinder->getParent());
        f->retractPlusEnd();
    }
};

/// Callback to retract the minus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionMinusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentRetractionMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = (Filament*)(_cylinder->getParent());
        f->retractMinusEnd();
    }
};

/// Callback to polymerize the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentPolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizePlusEnd();
    }
};

/// Callback to polymerize the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentPolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizeMinusEnd();
    }
};

/// Callback to depolymerize the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentDepolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = (Filament*)(_cylinder->getParent());
        
        f->depolymerizePlusEnd();
    }
};

/// Callback to depolymerize the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentDepolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = (Filament*)(_cylinder->getParent());
        f->depolymerizeMinusEnd();
    }
};

/// Callback to unbind a BranchingPoint from a Filament
struct BranchingPointUnbindingCallback {
    
    SubSystem* _ps;
    BranchingPoint* _branchingPoint;
    
    BranchingPointUnbindingCallback(BranchingPoint* b, SubSystem* ps)
    : _ps(ps), _branchingPoint(b) {}
    
    void operator() (ReactionBase *r) {
        
        //remove the branching point
        _ps->removeTrackable<BranchingPoint>(_branchingPoint);
        delete _branchingPoint;
    }
};


/// Callback to create a BranchingPoint on a Filament
struct BranchingCallback {
    
    SubSystem* _ps;        ///< ptr to subsystem
    
    BranchingManager* _bManager; ///< Branching manager for this compartment
    
    short _plusEnd;        ///< Plus end marker of new cylinder
    
    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction
    
    BranchingCallback(BranchingManager* bManager, short plusEnd,
                      float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _bManager(bManager),
    _plusEnd(plusEnd), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        BranchingPoint* b;
        float frate;
        short branchType = _bManager->getBoundInt();
        
        //choose a random binding site from manager
        auto site = _bManager->chooseBindingSite();
        
        //get info from site
        Cylinder* c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();
        
        double pos = double(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];
        if(SysParams::RUNSTATE==true){
        //Get a position and direction of a new filament
        auto x1 = c1->getFirstBead()->coordinate;
        auto x2 = c1->getSecondBead()->coordinate;
        
        //get original direction of cylinder
        auto p= midPointCoordinate(x1, x2, pos);
        vector<double> n = twoPointDirection(x1, x2);
        
        //get branch projection
#ifdef MECHANICS
        //use mechanical parameters
        double l, t;
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
        double l = 10.0;
        double t = 1.22;
#endif
        double s = SysParams::Geometry().monomerSize[filType];
        
        auto branchPosDir = branchProjection(n, p, l, s, t);
        auto bd = get<0>(branchPosDir); auto bp = get<1>(branchPosDir);
        
        //create a new filament
        Filament* f = _ps->addTrackable<Filament>(_ps, filType, bp, bd, true, true);
        
        //mark first cylinder
        Cylinder* c = f->getCylinderVector().front();
        c->getCCylinder()->getCMonomer(0)->speciesPlusEnd(_plusEnd)->up();
        
        //create new branch
            CMonomer* x=c->getCCylinder()->getCMonomer(0);
            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().brancherBoundIndex[filType]);
//                auto yy =c->getCCylinder()->getCMonomer(p)->speciesBrancher(branchType);
//                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                            }
            std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
            
        b= _ps->addTrackable<BranchingPoint>(c1, c, branchType, pos);
            
            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().brancherBoundIndex[filType]);
//                auto yy =c->getCCylinder()->getCMonomer(p)->speciesBrancher(branchType);
//                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
            }
            //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
        frate=_offRate;
        }
        else
        {
            CCylinder* c; auto check = false;
        vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> BrT=_bManager->getbtuple();
            for(auto T:BrT){
                CCylinder* cx=get<0>(get<0>(T));
                double p = double(get<1>(get<0>(T)))/ double(SysParams::Geometry().cylinderNumMon[filType]);
                if(cx->getCylinder()->getID()==c1->getID() && p==pos){
                    c=get<0>(get<1>(T));
                    check = true;
                    break;
                }}
            if(check){
///                auto cyl = c->getCylinder();
                //std::cout<<twoPointDistance(cyl->getFirstBead()->coordinate,cyl->getSecondBead()->coordinate)<<" ";
            //std::cout<<c->getCylinder()->getID()<<endl;
                CMonomer* x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().brancherBoundIndex[filType]);
//                    auto yy =c->getCMonomer(p)->speciesBrancher(branchType);
//                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
                
            b= _ps->addTrackable<BranchingPoint>(c1, c->getCylinder(), branchType, pos);
                
                x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().brancherBoundIndex[filType]);
//                    auto yy =c->getCMonomer(p)->speciesBrancher(branchType);
//                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
            frate=0.0;
            }
            else
                cout<<"Brancher Error. Cannot find binding Site in the list. Cannot complete restart. Exiting." <<endl;
                //exit(EXIT_FAILURE);
        }
        
        //create off reaction
        auto cBrancher = b->getCBranchingPoint();
        
        cBrancher->setRates(_onRate, frate);
        cBrancher->createOffReaction(r, _ps);
        cBrancher->getOffReaction()->setBareRate(SysParams::BUBBareRate[branchType]);
    }
};

/// Callback to unbind a CaMKIIingPoint from a Filament
struct CaMKIIingPointUnbindingCallback {
    
    SubSystem* _ps;
    CaMKIIingPoint* _camkiiingPoint;

    CaMKIIingPointUnbindingCallback(CaMKIIingPoint* b, SubSystem* ps)
    : _ps(ps), _camkiiingPoint(b) {}
    
    void operator() (ReactionBase *r) {
        cout<<"========== CaMKII Unbinding CallBack "; //Carlos verbose prints
        cout<< "ID: " << _camkiiingPoint->getID() << " Coord:" << _camkiiingPoint->getCoordinationNumber() << endl; //Carlos verbose prints

        // Removing the internal reaction of the CaMKII cylinder
        _camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);

        //remove the camkiiing point
        _camkiiingPoint->getCCaMKIIingPoint()->removeBond(_camkiiingPoint->getBond(0));
        cerr << "====== MILLAD: Removing CaMKIICylinder " << _camkiiingPoint->getCaMKIICylinder()  << "   " << _camkiiingPoint->getCaMKIICylinder()->getCCylinder() << endl;
        _ps->removeTrackable<CaMKIICylinder>(_camkiiingPoint->getCaMKIICylinder());
        _ps->removeTrackable<CaMKIIingPoint>(_camkiiingPoint);
        delete _camkiiingPoint;
    }
};

/// Callback to unbind a CaMKIIingPoint from a Filament
struct CaMKIIingPointUnbundlingCallback {

    SubSystem* _ps;
    CaMKIIingPoint* _camkiiingPoint;
    CaMKIIBundlingManager* _bManager;

    CaMKIIingPointUnbundlingCallback(CaMKIIingPoint* b, SubSystem* ps, CaMKIIBundlingManager* bManager)
            : _ps(ps), _camkiiingPoint(b), _bManager(bManager) {}

    void operator() (ReactionBase *r) {
        cout<<"========== CaMKII Unbundling CallBack "; //Carlos verbose prints
        cout<< "ID: " << _camkiiingPoint->getID() << " Coord:" << _camkiiingPoint->getCoordinationNumber()<< "-->"; //Carlos verbose prints
        auto bond = _camkiiingPoint->removeRandomBond();
        // Adding the cylinder assigned to the removed bond?
//        _bManager->addPossibleBindings(get<0>(bond)->getCCylinder(), get<1>(bond));
        cout<< _camkiiingPoint->getCoordinationNumber() <<endl;
        if(_camkiiingPoint->getCoordinationNumber() == 1L) {

            // passivate unbundling reaction
			r->passivateReaction();

			// Getting a handler to the current off-reaction
            ReactionBase *offRxnBinding = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();

            // Creating an off-reaction if the off-reaction does not exist
            // due to CaMKII moving between compartments
			if(offRxnBinding == nullptr) {
			    _camkiiingPoint->getCCaMKIIingPoint()->createOffReactionBinding(nullptr, _ps);
                ReactionBase *offRxn = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();
                _camkiiingPoint->getCCaMKIIingPoint()->setOffReaction(offRxn);
                offRxnBinding = offRxn;
			}

			// Reactivating the unbinding reaction
			offRxnBinding->activateReaction();

			// Removing the current off-reaction from the internal list in CaMKII cylinder
//			get<0>(_camkiiingPoint->getBonds()[0])->getCCylinder()->removeInternalReaction(r);
			_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);
//            get<0>(_camkiiingPoint->getBonds()[0])->getCCylinder()->addInternalReaction(offRxnBinding);
			_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->addInternalReaction(offRxnBinding);
            _camkiiingPoint->getCCaMKIIingPoint()->setOffReaction(offRxnBinding);

            // Removing the unbundling reaction
            _camkiiingPoint->getCCaMKIIingPoint()->setOffRxnBundling(nullptr);

        } else if(_camkiiingPoint->getCoordinationNumber() == _bManager->getMaxCoordinationNumber()-1 ) {
            _bManager->addPossibleBindings(_camkiiingPoint->getCaMKIICylinder()->getCCylinder());
        }

        _camkiiingPoint->updateReactionRates();
    }
};


/// Callback to create a CaMKIIingPoint on a Filament
struct CaMKIIBindingCallback {

    SubSystem* _ps;        ///< ptr to subsystem

    CaMKIIBindingManager* _bManager; ///< CaMKIIing manager for this compartment

    //short _plusEnd;        ///< Plus end marker of new cylinder

    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction

    CaMKIIBindingCallback(CaMKIIBindingManager* bManager, float onRate, float offRate, SubSystem* ps)
    : _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}

    void operator() (ReactionBase *r) {
        CaMKIIingPoint *b;
        cout << "========== CaMKII Binding CallBack "; //Carlos verbose prints
        float frate;
        short camkiiType = _bManager->getBoundInt();

        //choose a random binding site from manager
        auto site = _bManager->chooseBindingSite();
        //get info from site
        Cylinder *c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();

        double pos = double(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];

        //Get a position and direction of a new filament
        auto x1 = c1->getFirstBead()->coordinate;
        auto x2 = c1->getSecondBead()->coordinate;

        //get original direction of cylinder
        auto p = midPointCoordinate(x1, x2, pos);
        vector<double> n = twoPointDirection(x1, x2);

        //get camkii projection
#ifdef MECHANICS
        //use mechanical parameters
        double l, t;
        if (SysParams::Mechanics().CaMKIIStretchingL.size() != 0) {
            l = SysParams::Mechanics().CaMKIIStretchingL[camkiiType];
//            t = SysParams::Mechanics().CaMKIIBendingTheta[camkiiType];
            t = M_PI_2;
        } else {
//            cout << "CaMKIIing initialization cannot occur unless mechanical parameters are specified."
//            << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
//            << endl;
            l = 10.0;
            t = M_PI_2;
        }
#else
        cout << "CaMKIIing initialization cannot occur unless mechanics is enabled. Using"
			<< " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
			<< endl;
			double l = 10.0;
			double t = 1.22;
#endif
        double s = SysParams::Geometry().monomerSize[filType];

        // initial position of CaMKIIingPoint
        auto cp = camkiiProjection(n, p, l);

        if (SysParams::RUNSTATE == true) {

            b = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);
            cout << "ID: " << b->getID() << endl; //Carlos verbose prints
            frate = _offRate;
        } else {
            CCylinder *c;
            auto check = false;
            vector<tuple<tuple<CCylinder *, short>, tuple<CCylinder *, short>>> BrT = _bManager->getbtuple();
            for (auto T:BrT) {
                CCylinder *cx = get<0>(get<0>(T));
                double p = double(get<1>(get<0>(T))) / double(SysParams::Geometry().cylinderNumMon[filType]);
                if (cx->getCylinder()->getID() == c1->getID() && p == pos) {
                    c = get<0>(get<1>(T));
                    check = true;
                    break;
                }
            }
            if (check) {
                CMonomer *x = c->getCMonomer(0);
                vector<Cylinder *> cy{c1, c->getCylinder()};
                b = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);

                x = c->getCMonomer(0);
                frate = 0.0;
            } else {
                cout << "CaMKIIer Error. Cannot find binding Site in the list. Cannot complete restart. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
        }

        //create off reaction
        auto cCaMKIIer = b->getCCaMKIIingPoint();

        cCaMKIIer->setRates(_onRate, frate);
        cCaMKIIer->createOffReactionCaMKII(r, _ps, nullptr);
        cCaMKIIer->getOffReaction()->setBareRate(SysParams::CaMKIIUnbindingBareRate[camkiiType]);
    }
};

/// Callback to create a CaMKIIingPoint on a Filament
struct CaMKIIBundlingCallback {
    
    SubSystem* _ps;        ///< ptr to subsystem
    
    CaMKIIBundlingManager* _bManager; ///< CaMKIIing manager for this compartment
    
    //short _plusEnd;        ///< Plus end marker of new cylinder
    
    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction
    
    CaMKIIBundlingCallback(CaMKIIBundlingManager* bManager, float onRate, float offRate, SubSystem* ps)
    //TODO add minSearchDist maxSearchDist maxcoordNumber as local data members
    : _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
    	  CaMKIIingPoint* cp;
    	// TODO CAMKII either find cp from NN or to get the one the reaction was on it
//        float frate;
        cout<<"========== CaMKII Bundling CallBack ";
        short camkiiType = _bManager->getBoundInt();

        //choose a random binding site from manager
        //TODO Find a CAMKII with coordination number less than max number
        auto site = _bManager->chooseBindingSites();

        Cylinder*       c1  = get<0>(get<0>(site))->getCylinder();
        short           pos1 = get<1>(get<0>(site));

        Cylinder*       c2  = get<0>(get<1>(site))->getCylinder();
        short           pos2 = get<1>(get<1>(site));

        bool b1 = (dynamic_cast<CaMKIICylinder*>(c1)!= nullptr);
        bool b2 = (dynamic_cast<CaMKIICylinder*>(c2)!= nullptr);

        // If none of them or both of them are CaMKIIPoint, do nothing!
        if((!b1 && !b2) || (b1 && b2))
        	return;

        if(b1) {
            cp = dynamic_cast<CaMKIICylinder *>(c1)->getCaMKIIPointParent();

            // cp should be CaMKIIPoint
            cp->addBond(c2, pos2);
            cout<< "(b1) ID: " << cp->getID() << " Coord:" << cp->getCoordinationNumber()<< " CaMKIItype:" << camkiiType <<endl; //Carlos verbose prints
        } else {
            cp = dynamic_cast<CaMKIICylinder *>(c2)->getCaMKIIPointParent();

            // cp should be CaMKIIPoint
            cp->addBond(c1, pos1);
            cout<< "(b2) ID: " << cp->getID() << " Coord:" << cp->getCoordinationNumber()<< " CaMKIItype:" << camkiiType <<endl; //Carlos verbose prints
        }

//TODO CJY make sure isn't needed before cleaning
#if 0

        // CamKIIingPoint cylinder pair filament with

        //get info from site
        Cylinder* c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();

        double pos = double(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];
        if(SysParams::RUNSTATE==true){
        //Get a position and direction of a new filament
        auto x1 = c1->getFirstBead()->coordinate;
        auto x2 = c1->getSecondBead()->coordinate;

        //get original direction of cylinder
        auto p= midPointCoordinate(x1, x2, pos);
        vector<double> n = twoPointDirection(x1, x2);

        //get camkii projection
#ifdef MECHANICS
        //use mechanical parameters
        double l, t;
        if(SysParams::Mechanics().CaMKIIStretchingL.size() != 0) {
            l = SysParams::Mechanics().CaMKIIStretchingL[camkiiType];
            t = SysParams::Mechanics().CaMKIIBendingTheta[camkiiType];
        }
        else {
            // TODO
            cout << "CaMKIIing  bundling initialization cannot occur unless mechanical parameters are specified."
            << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
            << endl;
            l = 10.0;
            t = 1.22;
        }
#else
        cout << "CaMKIIing initialization cannot occur unless mechanics is enabled. Using"
        << " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
        << endl;
        double l = 10.0;
        double t = 1.22;
#endif
        double s = SysParams::Geometry().monomerSize[filType];

        auto camkiiPosDir = camkiiProjection(n, p, l, s, t);
        auto bd = get<0>(camkiiPosDir); auto bp = get<1>(camkiiPosDir);

        //TODO Search for new filament
//        //create a new filament
//        Filament* f = _ps->addTrackable<Filament>(_ps, filType, bp, bd, true, true);
//
//        //mark first cylinder

//        auto bindingSite = _bManager->chooseBindingSite();
//        auto ccx = get<0>(bindingSite);
//        auto cylin = ccx->getCylinder();
//        short bs = get<1>(ccx);

//        Cylinder* c = f->getCylinderVector().front();
//        //c->getCCylinder()->getCMonomer(0)->speciesPlusEnd(_plusEnd)->up();
//
//        //create new camkii
//            CMonomer* x=c->getCCylinder()->getCMonomer(0);
//            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
//                auto yy =c->getCCylinder()->getCMonomer(p)->speciesCaMKIIer(camkiiType);
//                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
//                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
//                            }
//            std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
//            vector<Cylinder*> cy{c1,c};
        //TODO Update new CAMKII to update coordination number
        //TODO Create the bonds between CAMKII and filament
        //b= _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos);

            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                auto yy =c->getCCylinder()->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
            }
            //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
        frate=_offRate;
        }
        else
        {
            CCylinder* c; auto check = false;
        vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> BrT=_bManager->getbtuple();
            for(auto T:BrT){
                CCylinder* cx=get<0>(get<0>(T));
                double p = double(get<1>(get<0>(T)))/ double(SysParams::Geometry().cylinderNumMon[filType]);
                if(cx->getCylinder()->getID()==c1->getID() && p==pos){
                    c=get<0>(get<1>(T));
                    check = true;
                    break;
                }}
            if(check){
                auto cyl = c->getCylinder();
                //std::cout<<twoPointDistance(cyl->getFirstBead()->coordinate,cyl->getSecondBead()->coordinate)<<" ";
            //std::cout<<c->getCylinder()->getID()<<endl;
                CMonomer* x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                    auto yy =c->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
                //vector<Cylinder*> cy{c1,c->getCylinder()};
                b= _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos);

                x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                    auto yy =c->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
            frate=0.0;
            }
            else
                cout<<"CaMKIIer Error. Cannot find binding Site in the list. Cannot complete restart. Exiting." <<endl;
                //exit(EXIT_FAILURE);
        }
#endif

        //create off reaction
        auto cCaMKIIer = cp->getCCaMKIIingPoint();
        cCaMKIIer->setRates(_onRate, _offRate);
        cCaMKIIer->createOffReactionCaMKII(r, _ps, _bManager);
        cCaMKIIer->getOffReaction()->setBareRate(SysParams::CaMKIIUnbundlingBareRate[camkiiType]);
        _bManager->updateCaMKIIPossibleBindings(cp->getCaMKIICylinder()->getCCylinder());
    }
};

/// Callback to unbind a Linker from a Filament
struct LinkerUnbindingCallback {
    
    SubSystem* _ps;
    Linker* _linker;
    
    LinkerUnbindingCallback(Linker* l, SubSystem* ps) : _ps(ps), _linker(l) {}
    
    void operator() (ReactionBase *r) {
        
        //remove the linker
        _ps->removeTrackable<Linker>(_linker);
        delete _linker;
    }
};

/// Callback to bind a Linker to Filament
struct LinkerBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    LinkerBindingManager* _lManager; ///< Linker binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction
    
    LinkerBindingCallback(LinkerBindingManager* lManager,
                          float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _lManager(lManager), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        //get a random binding
        short linkerType = _lManager->getBoundInt();
        float f;
        
        //choose a random binding site from manager
        auto site = _lManager->chooseBindingSites();
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        short filType = c1->getType();
        
        // Create a linker
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        
        double pos1 = double(get<1>(site[0])) / cylinderSize;
        double pos2 = double(get<1>(site[1])) / cylinderSize;
        
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
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        //aravind june 24, 2016
        cLinker->getOffReaction()->setBareRate(SysParams::LUBBareRate[linkerType]);
        //@
        l->updateReactionRates();
#endif
    }
};

/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    MotorGhost* _motor;
    
    MotorUnbindingCallback(MotorGhost* m, SubSystem* ps) :
    
    _ps(ps), _motor(m) {}
    
    void operator() (ReactionBase *r) {
        
        //find the motor binding manager associated with this species
//        Species* sd = &(r->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
//        Compartment* c = static_cast<Compartment*>(sd->getParent());
//        auto mManager = c->getMotorBindingManager(_motor->getType());
        
        //DEPRECATED AS OF 9/8/16
//        
//        mManager->removeUnboundID(MotorGhost::_motorGhosts.deleteID());
//        
//        //re-add unbound ID to motor binding manager
//        mManager->addUnboundID(_motor->getID());
        
        //remove the motor
        _ps->removeTrackable<MotorGhost>(_motor);
        delete _motor;
    }
};


/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    MotorBindingManager* _mManager;///< Motor binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction
    
    MotorBindingCallback(MotorBindingManager* mManager,
                         float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _mManager(mManager), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        //get a random binding
        short motorType = _mManager->getBoundInt();
        float f;
        //choose a random binding site from manager
        auto site = _mManager->chooseBindingSites();
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        short filType = c1->getType();
        
        // Create a motor
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        
        double pos1 = double(get<1>(site[0])) / cylinderSize;
        double pos2 = double(get<1>(site[1])) / cylinderSize;
        
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
    
    MotorWalkingCallback(Cylinder* c,
                         short oldPosition, short newPosition,
                         short motorType, short boundType, SubSystem* ps)
    
    :_c(c), _oldPosition(oldPosition), _newPosition(newPosition),
    _motorType(motorType), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //get species
        CCylinder* cc = _c->getCCylinder();
        CMonomer* monomer = cc->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        
        short filType = _c->getType();
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        double oldpos = double(_oldPosition) / cylinderSize;
        double newpos = double(_newPosition) / cylinderSize;
        
        m->moveMotorHead(_c, oldpos, newpos, _boundType, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        m->updateReactionRates();

#endif
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
    
    MotorMovingCylinderCallback(Cylinder* oldC, Cylinder* newC,
                                short oldPosition, short newPosition,
                                short motorType, short boundType, SubSystem* ps)
    
    :_oldC(oldC), _newC(newC), _oldPosition(oldPosition), _newPosition(newPosition),
    _motorType(motorType), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //get species
        CCylinder* oldCC = _oldC->getCCylinder();
        CMonomer* monomer = oldCC->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        short filType = _oldC->getType();
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
        double oldpos = double(_oldPosition) / cylinderSize;
        double newpos = double(_newPosition) / cylinderSize;
        
        m->moveMotorHead(_oldC, _newC, oldpos, newpos, _boundType, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        m->updateReactionRates();
#endif
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
    
    FilamentCreationCallback(short plusEnd, short minusEnd, short filament,
                             short filType, SubSystem* ps, Compartment* c = nullptr)
    
    : _plusEnd(plusEnd), _minusEnd(minusEnd), _filament(filament),
    _filType(filType), _compartment(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        Compartment* c;
        
        //no compartment was set, pick a random one
        if(_compartment == nullptr)
            c = GController::getRandomCompartment();
        else
            c = _compartment;
        
        //set up a random initial position and direction
        vector<double> position;
        vector<double> direction;
        
        while(true) {
            position = GController::getRandomCoordinates(c);
            
            //getting random numbers between -1 and 1
            direction = {Rand::randDouble(-1,1), Rand::randDouble(-1,1), Rand::randDouble(-1,1)};
            normalize(direction);
            
            auto npp = nextPointProjection(position,
                                           SysParams::Geometry().cylinderSize[_filType], direction);
            
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
    
};

///Struct to sever a filament based on a reaction
struct FilamentSeveringCallback {
    
    Cylinder* _c1;  ///< Filament severing point
    
    FilamentSeveringCallback(Cylinder* c1) : _c1(c1) {}
    
    void operator() (ReactionBase* r) {
        
        //reactants should be re-marked
        for(int i = 0; i < SEVERINGREACTANTS + 1; i++)
            r->rspecies()[i]->up();
        
        //sever the filament at given position
        Filament* f = (Filament*)(_c1->getParent());
        
        //create a new filament by severing
        Filament* newF = f->sever(_c1->getPosition());
        
        //if severing did not occur, return
        if(newF == nullptr) return;
    }
};

/// Struct to destroy a filament based on a reaction
struct FilamentDestructionCallback {
    
    Cylinder* _c; ///< Cylinder to destroy
    
    SubSystem* _ps; ///< SubSystem ptr
    
    FilamentDestructionCallback(Cylinder* c, SubSystem* ps) : _c(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        Filament* f = (Filament*)(_c->getParent());
        
        _ps->removeTrackable<Filament>(f);
        delete f;
    }
};

#endif

#endif
