
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

#ifndef M3SYM_ChemCallbacks_h
#define M3SYM_ChemCallbacks_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Boundary.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

#ifdef REACTION_SIGNALING

/// Callback to extend the front of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionFrontCallback {
    
    Cylinder* _cylinder;
    
    short _plusEnd; ///< Plus end species to mark
    
    //Constructor, sets members
    FilamentExtensionFrontCallback(Cylinder* cylinder, short plusEnd)
        : _cylinder(cylinder), _plusEnd(plusEnd){};
    
    //Callback
    void operator() (ReactionBase *r){
        //extend the front
        Filament* f = _cylinder->getFilament();
        f->extendFront(_plusEnd);
    }
};

/// Callback to extend the back of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionBackCallback {
    
    Cylinder* _cylinder;
    
    short _minusEnd; ///< Minus end species to mark
    
    //Constructor, sets members
    FilamentExtensionBackCallback(Cylinder* cylinder, short minusEnd)
        : _cylinder(cylinder), _minusEnd(minusEnd){};
    //Callback
    void operator() (ReactionBase *r){
        //extend the back
        Filament* f = _cylinder->getFilament();
        f->extendBack(_minusEnd);
    }
};

/// Callback to retract the front of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionFrontCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentRetractionFrontCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->retractFront();
    }
};

/// Callback to retract the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionBackCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentRetractionBackCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->retractBack();
    }
};

/// Callback to polymerize the front of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationFrontCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentPolymerizationFrontCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->polymerizeFront();
    }
};

/// Callback to polymerize the back of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationBackCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentPolymerizationBackCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->polymerizeBack();
    }
};

/// Callback to depolymerize the front of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationFrontCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentDepolymerizationFrontCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->depolymerizeFront();
    }
};

/// Callback to depolymerize the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationBackCallback {
    
    Cylinder* _cylinder;
    
    //Constructor, sets members
    FilamentDepolymerizationBackCallback(Cylinder* cylinder)
        : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
        Filament* f = _cylinder->getFilament();
        f->depolymerizeBack();
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
        _ps->removeBranchingPoint(_branchingPoint);
    }
};

/// Callback to create a BranchingPoint on a Filament
struct BranchingPointCreationCallback {
    
    SubSystem* _ps;        ///< ptr to subsystem
    Cylinder* _c1;         ///< Cylinders to attach this branchingpoint to
    
    short _branchType;     ///< Type of branchingpoint
    short _position;       ///< Position to attach this branchingpoint
    
    short _plusEnd;        ///< Plus end marker of new cylinder
    
    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction
    
    BranchingPointCreationCallback(Cylinder* c1,
                                   short branchType,
                                   short plusEnd,
                                   short position,
                                   float onRate,
                                   float offRate,
                                   SubSystem* ps)
    
        : _ps(ps), _c1(c1), _branchType(branchType), _position(position),
          _plusEnd(plusEnd), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        double pos = double(_position) / SysParams::Geometry().cylinderIntSize;
        
        //Get a position and direction of a new filament
        auto x1 = _c1->getFirstBead()->coordinate;
        auto x2= _c1->getSecondBead()->coordinate;
        
        //get original direction of cylinder
        auto p= midPointCoordinate(x1, x2, pos);
        vector<double> n = twoPointDirection(x1, x2);
        
        //get branch projection
#ifdef MECHANICS
        //use mechanical parameters
        double l = SysParams::Mechanics().BrStretchingL[_branchType];
        double t = SysParams::Mechanics().BrBendingTheta[_branchType];
#else
        cout << "Branching reaction cannot occur unless mechanics is enabled. Using"
             << " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
             << endl;
        double l = 10.0;
        double t = 1.22;
#endif
        double s = SysParams::Geometry().monomerSize;
        
        auto branchPosDir = branchProjection(n, p, l, s, t);
        auto bd = get<0>(branchPosDir);
        auto bp = get<1>(branchPosDir);
        
        //create a new filament
        Filament* f = _ps->addNewFilament(bp, bd, true);
        
        //mark first cylinder
        Cylinder* c = f->getCylinderVector().front();
        CMonomer* m = c->getCCylinder()->getCMonomer(0);
        m->speciesPlusEnd(_plusEnd)->up();
        
        //create new branch
        BranchingPoint* b= _ps->addNewBranchingPoint(_c1, c, _branchType, pos);
        
        //create off reaction
        b->getCBranchingPoint()->setOnRate(_onRate);
        b->getCBranchingPoint()->setOffRate(_offRate);
        b->getCBranchingPoint()->createOffReaction(r, _ps);
    }
};

/// Callback to unbind a Linker from a Filament
struct LinkerUnbindingCallback {
    
    SubSystem* _ps;
    Linker* _linker;
    
    LinkerUnbindingCallback(Linker* l, SubSystem* ps) : _ps(ps), _linker(l) {}
    
    void operator() (ReactionBase *r) {
        //remove the linker
        _ps->removeLinker(_linker);
    }
};

/// Callback to bind a Linker to Filament
struct LinkerBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    Cylinder* _c1, *_c2;          ///< Cylinders to attach this linker to
    
    short _linkerType;            ///< Type of linker
    short _position1, _position2; ///< Positions to attach this linker
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    LinkerBindingCallback(Cylinder* c1, Cylinder* c2,
                          short linkerType,
                          short position1,
                          short position2,
                          float onRate,
                          float offRate,
                          SubSystem* ps)
    
        : _ps(ps), _c1(c1), _c2(c2), _linkerType(linkerType),
          _position1(position1), _position2(position2),
          _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        // Create a linker
        int cylinderSize = SysParams::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        Linker* l = _ps->addNewLinker(_c1, _c2, _linkerType, pos1, pos2);
        
        //create off reaction
        l->getCLinker()->setOnRate(_onRate);
        l->getCLinker()->setOffRate(_offRate);
        l->getCLinker()->createOffReaction(r, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        l->updateReactionRates();
#endif
    }
};

/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    MotorGhost* _motor;
    
    MotorUnbindingCallback(MotorGhost* m, SubSystem* ps) : _ps(ps), _motor(m) {}
    
    void operator() (ReactionBase *r) {
        //remove the motor
        _ps->removeMotorGhost(_motor);
    }
};


/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;               ///< Ptr to subsystem
    Cylinder* _c1, *_c2;          ///< Cylinders to attach this motor to
    
    short _motorType;             ///< Type of motor
    short _position1, _position2; ///< Positions to attach this motor
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction
    
    MotorBindingCallback(Cylinder* c1, Cylinder* c2,
                         short motorType,
                         short position1,
                         short position2,
                         float onRate,
                         float offRate,
                         SubSystem* ps)
    
        : _ps(ps), _c1(c1), _c2(c2), _motorType(motorType),
          _position1(position1), _position2(position2),
          _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {

        // Create a motor
        int cylinderSize = SysParams::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        MotorGhost* m = _ps->addNewMotorGhost(_c1, _c2, _motorType, pos1, pos2);

        //create off reaction
        m->getCMotorGhost()->setOnRate(_onRate);
        m->getCMotorGhost()->setOffRate(_offRate);
        m->getCMotorGhost()->createOffReaction(r, _ps);
        
#ifdef DYNAMICRATES
        //reset the associated walking reactions
        m->updateReactionRates();
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
                        short oldPosition,
                        short newPosition,
                        short motorType,
                        short boundType,
                        SubSystem* ps)
    
        :_c(c), _oldPosition(oldPosition), _newPosition(newPosition),
         _motorType(motorType), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //get species
        CCylinder* cc = _c->getCCylinder();
        CMonomer* monomer = cc->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderIntSize;
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
                                short oldPosition,
                                short newPosition,
                                short motorType,
                                short boundType,
                                SubSystem* ps)
    
        :_oldC(oldC), _newC(newC), _oldPosition(oldPosition), _newPosition(newPosition),
         _motorType(motorType), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //get species
        CCylinder* oldCC = _oldC->getCCylinder();
        CMonomer* monomer = oldCC->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(_motorType);
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderIntSize;
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

    Compartment* _compartment; ///< compartment to put this filament in
    SubSystem* _ps; ///< Ptr to the subsystem
    
    FilamentCreationCallback(short plusEnd,
                             short minusEnd,
                             short filament,
                             SubSystem* ps,
                             Compartment* c = nullptr)
    
        : _plusEnd(plusEnd), _minusEnd(minusEnd), _filament(filament),
          _compartment(c), _ps(ps) {}
    
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
            direction = {randomDouble(-1,1),
                         randomDouble(-1,1),
                         randomDouble(-1,1)};
            
            auto npp = nextPointProjection(position,
                SysParams::Geometry().cylinderSize, direction);
            
            //check if within boundary
            if(_ps->getBoundary()->within(position) &&
               _ps->getBoundary()->within(npp))
                break;
        }
        
        //create filament, set up ends and filament species
        Filament* f = _ps->addNewFilament(position, direction);
        
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
        Filament* f = _c1->getFilament();
        
        //create a new filament by severing
        f->sever(_c1->getPositionFilament());
    }
};

/// Struct to destroy a filament based on a reaction
struct FilamentDestructionCallback {
    
    Cylinder* _c; ///< Cylinder to destroy
    
    FilamentDestructionCallback(Cylinder* c) : _c(c) {}
    
    void operator() (ReactionBase* r) {
        Filament* f = _c->getFilament();
        delete f;
    }
};

#endif

#endif
