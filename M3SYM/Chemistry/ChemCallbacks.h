
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
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
#include "ChemRNode.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

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
        f->extendFront();
        
        //get last cylinder, mark species
        auto newCylinder = f->getCylinderVector().back();
        CMonomer* m = newCylinder->getCCylinder()->getCMonomer(0);
        
        m->speciesPlusEnd(_plusEnd)->up();
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
        f->extendBack();
        
        //get first cylinder, mark species
        auto newCylinder = f->getCylinderVector().front();
        auto newCCylinder = newCylinder->getCCylinder();
        CMonomer* m = newCCylinder->getCMonomer(newCCylinder->getSize() - 1);
        
        m->speciesMinusEnd(_minusEnd)->up();
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
        
#ifdef DYNAMICRATES
        //update rates of new front
        f->getCylinderVector().back()->updateReactionRates();
#endif
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
        
#ifdef DYNAMICRATES
        //update rates of new back
        f->getCylinderVector().front()->updateReactionRates();
#endif
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
        //remove the unbinding reaction
        CCylinder* cc = _branchingPoint->getFirstCylinder()->getCCylinder();
        cc->removeInternalReaction(r);
        
        //mark the correct species on the minus end of the branched
        //filament. If this is a filament species, change it to its
        //corresponding minus end. If a plus end, release a diffusing
        //or bulk species, depending on the initial reaction.
        CCylinder* childCC = _branchingPoint->getSecondCylinder()->getCCylinder();
        CMonomer* m = childCC->getCMonomer(0);
        short speciesFilament = m->activeSpeciesFilament();
        
        //there is a filament species, mark its corresponding minus end
        if(speciesFilament != -1) {
            m->speciesMinusEnd(speciesFilament)->up();
        }
        //mark the free species instead
        else {
            //find the free species
            short speciesPlusEndNum = m->activeSpeciesPlusEnd();
            Species* speciesFilament = m->speciesFilament(speciesPlusEndNum);
            
            string speciesName = SpeciesNamesDB::Instance()->
                                 removeUniqueName(speciesFilament->getName());
            Species* freeMonomer =
                _branchingPoint->getCompartment()->findSpeciesByName(speciesName);
            
            //remove the filament from the system
            delete _branchingPoint->getSecondCylinder()->getFilament();
            
            //update reaction rates
            for(auto &r : freeMonomer->getRSpecies().reactantReactions())
                r->getRNode()->activateReaction();
        }
        
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
    
    float _offRate;        ///< Rate of the unbinding reaction
    
    BranchingPointCreationCallback(Cylinder* c1, short branchType, short plusEnd,
                                   short position, float offRate, SubSystem* ps)
        : _ps(ps), _c1(c1), _branchType(branchType), _plusEnd(plusEnd),
          _position(position),  _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        double pos = double(_position) / cylinderSize;
        
        //Get a position and direction of a new filament
        auto position1 = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                            _c1->getSecondBead()->coordinate, pos);
        
        //randomize position and direction
        double rand1 = randomDouble(-1,1);
        double rand2 = randomDouble(-1,1);
        double rand3 = randomDouble(-1,1);
        position1[0] = position1[0] + rand1;
        position1[1] = position1[1] + rand2;
        position1[2] = position1[2] + rand3;
        
        auto position2 = _c1->getSecondBead()->coordinate;
        position2[0] = position2[0] + 2*rand1;
        position2[1] = position2[1] + 2*rand2;
        position2[2] = position2[2] + 2*rand3;
        
        vector<double> direction = twoPointDirection(position1,position2);
        
        //create a new filament
        Filament* f = _ps->addNewFilament(position1, direction, true);
        
        //mark first cylinder
        Cylinder* c = f->getCylinderVector()[0];
        CMonomer* m = c->getCCylinder()->getCMonomer(0);
        m->speciesPlusEnd(_plusEnd)->up();
        
        //create new branch
        BranchingPoint* b= _ps->addNewBranchingPoint(_c1, c, _branchType, pos);
        
        //add the unbinding reaction and callback
        //first, find the correct diffusing or bulk species
        Reaction<BRANCHINGREACTANTS, BRANCHINGPRODUCTS - 1>* branchReact =
            dynamic_cast<Reaction<BRANCHINGREACTANTS, BRANCHINGPRODUCTS - 1>*>(r);
        Species* freeBrancher = &(branchReact->rspecies()[0]->getSpecies());
        
        //create the reaction species
        m = _c1->getCCylinder()->getCMonomer(_position);
        vector<Species*> offSpecies =
            {m->speciesBrancher(_branchType), m->speciesBound(0), freeBrancher};
        
        ReactionBase* offRxn =
        new Reaction<BUNBINDINGREACTANTS,BUNBINDINGPRODUCTS>(offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::BRANCHUNBINDING);
        
        BranchingPointUnbindingCallback bcallback(b, _ps);
        boost::signals2::shared_connection_block
            rcb(offRxn->connect(bcallback,false));
        
        _c1->getCCylinder()->addInternalReaction(offRxn);
        b->getCBranchingPoint()->setOffReaction(offRxn);
        
    }
};

/// Callback to unbind a Linker from a Filament
struct LinkerUnbindingCallback {
    
    SubSystem* _ps;
    Linker* _linker;
    
    LinkerUnbindingCallback(Linker* l, SubSystem* ps) : _ps(ps), _linker(l) {}
    
    void operator() (ReactionBase *r) {
        //remove the unbinding reaction
        CCylinder* cc1 = _linker->getFirstCylinder()->getCCylinder();
        CCylinder* cc2 = _linker->getSecondCylinder()->getCCylinder();
        cc1->removeCrossCylinderReaction(cc2, r);

#ifdef DYNAMICRATES
        //reset the associated reactions
        _linker->getMLinker()->stretchForce = 0.0;
        _linker->updateReactionRates();
#endif
        
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
    
    float _offRate;               ///< Rate of the unbinding reaction

    LinkerBindingCallback(Cylinder* c1, Cylinder* c2,
                          short linkerType,
                          short position1, short position2,
                          float offRate, SubSystem* ps)
        : _ps(ps), _c1(c1), _c2(c2), _linkerType(linkerType),
          _position1(position1), _position2(position2), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        // Create a linker
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        Linker* l = _ps->addNewLinker(_c1, _c2, _linkerType, pos1, pos2);
        
        //Attach the callback to the off reaction, add it
        //first, get species from this rxn
        Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>* onRxn =
            (Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>*)r;
        RSpecies** rSpecies = onRxn->rspecies();
        vector<Species*> offSpecies;
        
        //copy into offspecies vector in opposite order
        for(int i = LMBINDINGREACTANTS;
            i < LMBINDINGREACTANTS + LMBINDINGPRODUCTS; i++)
            offSpecies.push_back(&rSpecies[i]->getSpecies());
        for(int i = 0; i < LMBINDINGREACTANTS; i++)
            offSpecies.push_back(&rSpecies[i]->getSpecies());
        
        ReactionBase* offRxn =
        new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::LINKERUNBINDING);
        
        LinkerUnbindingCallback lcallback(l, _ps);
        boost::signals2::shared_connection_block
            rcb(offRxn->connect(lcallback,false));
        
        _c1->getCCylinder()->addCrossCylinderReaction(_c2->getCCylinder(), offRxn);
        l->getCLinker()->setOffReaction(offRxn);
    }
};

/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    MotorGhost* _motor;
    
    MotorUnbindingCallback(MotorGhost* m, SubSystem* ps) : _ps(ps), _motor(m) {}
    
    void operator() (ReactionBase *r) {
        //remove the unbinding reaction
        CCylinder* cc1 = _motor->getFirstCylinder()->getCCylinder();
        CCylinder* cc2 = _motor->getSecondCylinder()->getCCylinder();
        cc1->removeCrossCylinderReaction(cc2, r);
        
#ifdef DYNAMICRATES
        //reset the associated reactions
        _motor->getMMotorGhost()->stretchForce = 0.0;
        _motor->updateReactionRates();
#endif
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
    
    float _offRate;               ///< Rate of the unbinding reaction
    
    MotorBindingCallback(Cylinder* c1, Cylinder* c2,
                         short motorType,
                         short position1, short position2,
                         float offRate, SubSystem* ps)
        : _ps(ps), _c1(c1), _c2(c2),
          _motorType(motorType),
          _position1(position1), _position2(position2), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {

        // Create a motor
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        MotorGhost* m = _ps->addNewMotorGhost(_c1, _c2, _motorType, pos1, pos2);
        
        //Attach the callback to the off reaction, add it
        //first, get species from this rxn
        Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>* onRxn =
        (Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>*)r;
        RSpecies** rSpecies = onRxn->rspecies();
        vector<Species*> offSpecies;
        
        //copy into offspecies vector in opposite order
        for(int i = LMBINDINGREACTANTS;
            i < LMBINDINGREACTANTS + LMBINDINGPRODUCTS; i++)
            offSpecies.push_back(&rSpecies[i]->getSpecies());
        for(int i = 0; i < LMBINDINGREACTANTS; i++)
            offSpecies.push_back(&rSpecies[i]->getSpecies());
        
        ReactionBase* offRxn =
        new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
            rcb(offRxn->connect(mcallback,false));
        
        _c1->getCCylinder()->addCrossCylinderReaction(_c2->getCCylinder(), offRxn);
        m->getCMotorGhost()->setOffReaction(offRxn);
    }
};

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingForwardCallback {
    
    Cylinder* _c;        ///< Cylinder this callback is attached to
    
    short _oldPosition;  ///< Old position of motor head
    short _newPosition;  ///< New position of motor head
    
    short _motorType;    ///< Type of motor
    short _boundType;    ///< Type of bound this motor took place of
    
    SubSystem* _ps;      ///< Ptr to subsystem
    
    MotorWalkingForwardCallback(Cylinder* c,
                                short oldPosition, short newPosition,
                                short motorType, short boundType,
                                SubSystem* ps)
        :_c(c), _motorType(motorType), _boundType(boundType),
         _oldPosition(oldPosition), _newPosition(newPosition), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //Find the motor
        CCylinder* cc = _c->getCCylinder();
        
        SpeciesMotor* sm1 = cc->getCMonomer(_oldPosition)->speciesMotor(_motorType);
        SpeciesMotor* sm2 = cc->getCMonomer(_newPosition)->speciesMotor(_motorType);
        SpeciesBound* sb2 = cc->getCMonomer(_newPosition)->speciesBound(_boundType);
        
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double shift =  1.0 / SystemParameters::Chemistry().numBindingSites;
        double newPosition;
        ReactionBase* newOffRxn, *offRxn;
        
        if(m->getCMotorGhost()->getFirstSpecies() == sm1) {

            newPosition = m->getFirstPosition() + shift;
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
                     m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                        ({sm2, s1, sb2, s3, s4}, offRxn->getRate());
        }
        else {
            newPosition = m->getSecondPosition() + shift;
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
                     m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                        ({s0, sm2, s2, sb2, s4}, offRxn->getRate());
        }
        //set new reaction type
        newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        //attach signal
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
            rcb(newOffRxn->connect(mcallback,false));
        
        CCylinder* cc1 = m->getFirstCylinder()->getCCylinder();
        CCylinder* cc2 = m->getSecondCylinder()->getCCylinder();
        
        //remove old reaction, add new one
        cc1->removeCrossCylinderReaction(cc2, offRxn);
        cc1->addCrossCylinderReaction(cc2, newOffRxn);
    
        //set new unbinding reaction
        m->getCMotorGhost()->setOffReaction(newOffRxn);
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderForwardCallback {
    
    Cylinder* _oldC;        ///< Old cylinder the motor is attached to
    Cylinder* _newC;        ///< New cylinder motor will be attached to
    
    short _oldPosition;     ///< Old position of motor head
    short _newPosition;     ///< New position of motor head
    
    short _motorType;       ///< Type of motor
    short _boundType;       ///< Type of bound this motor is taking place of
    
    SubSystem* _ps;         ///< Ptr to subsystem
    
    MotorMovingCylinderForwardCallback(Cylinder* oldC, Cylinder* newC,
                                       short oldPosition,
                                       short newPosition,
                                       short motorType,
                                       short boundType,
                                       SubSystem* ps)
        :_oldC(oldC), _newC(newC), _motorType(motorType), _boundType(boundType),
         _oldPosition(oldPosition), _newPosition(newPosition), _ps(ps){}
    
    void operator() (ReactionBase* r) {
        
        //Find the motor
        CCylinder* oldCC = _oldC->getCCylinder();
        CCylinder* newCC = _newC->getCCylinder();
        
        SpeciesMotor* sm1 =oldCC->getCMonomer(_oldPosition)->speciesMotor(_motorType);
        SpeciesMotor* sm2 =newCC->getCMonomer(_newPosition)->speciesMotor(_motorType);
        SpeciesBound* sb2 =newCC->getCMonomer(_newPosition)->speciesBound(_boundType);
        
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        ReactionBase* newOffRxn, *offRxn;
        CCylinder* cc1, *cc2;
        //initial set of cylinders
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        
        if(m->getCMotorGhost()->getFirstSpecies() == sm1) {
            
            m->setFirstCylinder(_newC);
            m->setFirstPosition((double)_newPosition
               / SystemParameters::Geometry().cylinderIntSize);
            m->getCMotorGhost()->setFirstSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
                     m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                        ({sm2, s1, sb2, s3, s4}, offRxn->getRate());
            
            //remove old off reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        else {
            m->setSecondCylinder(_newC);
            m->setSecondPosition((double)_newPosition
               / SystemParameters::Geometry().cylinderIntSize);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
                     m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                        ({s0, sm2, s2, sb2, s4}, offRxn->getRate());
            
            //remove old off reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        
        //set new reaction type
        newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        //attach signal
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
            rcb(newOffRxn->connect(mcallback,false));
        
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        
        //add new
        cc1->addCrossCylinderReaction(cc2, newOffRxn);
        
        //set new unbinding reaction
        m->getCMotorGhost()->setOffReaction(newOffRxn);
        
    }
};

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingBackwardCallback {
    
    Cylinder* _c;        ///< Cylinder this callback is attached to
    
    short _oldPosition;  ///< Old position of motor head
    short _newPosition;  ///< New position of motor head
    
    short _motorType;    ///< Type of motor
    short _boundType;    ///< Type of bound this motor took place of
    
    SubSystem* _ps;      ///< Ptr to subsystem
    
    MotorWalkingBackwardCallback(Cylinder* c,
                                short oldPosition, short newPosition,
                                short motorType, short boundType,
                                SubSystem* ps)
    :_c(c), _motorType(motorType), _boundType(boundType),
    _oldPosition(oldPosition), _newPosition(newPosition), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        
        //Find the motor
        CCylinder* cc = _c->getCCylinder();
        
        SpeciesMotor* sm1 = cc->getCMonomer(_oldPosition)->speciesMotor(_motorType);
        SpeciesMotor* sm2 = cc->getCMonomer(_newPosition)->speciesMotor(_motorType);
        SpeciesBound* sb2 = cc->getCMonomer(_newPosition)->speciesBound(_boundType);
        
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double shift =  1.0 / SystemParameters::Chemistry().numBindingSites;
        double newPosition;
        ReactionBase* newOffRxn, *offRxn;
        
        if(m->getCMotorGhost()->getFirstSpecies() == sm1) {
            
            newPosition = m->getFirstPosition() - shift;
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
            m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
            ({sm2, s1, sb2, s3, s4}, offRxn->getRate());
        }
        else {
            newPosition = m->getSecondPosition() - shift;
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
            m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
            ({s0, sm2, s2, sb2, s4}, offRxn->getRate());
        }
        //set new reaction type
        newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        //attach signal
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
        rcb(newOffRxn->connect(mcallback,false));
        
        CCylinder* cc1 = m->getFirstCylinder()->getCCylinder();
        CCylinder* cc2 = m->getSecondCylinder()->getCCylinder();
        
        //remove old reaction, add new one
        cc1->removeCrossCylinderReaction(cc2, offRxn);
        cc1->addCrossCylinderReaction(cc2, newOffRxn);
        
        //set new unbinding reaction
        m->getCMotorGhost()->setOffReaction(newOffRxn);
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderBackwardCallback {
    
    Cylinder* _oldC;        ///< Old cylinder the motor is attached to
    Cylinder* _newC;        ///< New cylinder motor will be attached to
    
    short _oldPosition;     ///< Old position of motor head
    short _newPosition;     ///< New position of motor head
    
    short _motorType;       ///< Type of motor
    short _boundType;       ///< Type of bound this motor is taking place of
    
    SubSystem* _ps;         ///< Ptr to subsystem
    
    MotorMovingCylinderBackwardCallback(Cylinder* oldC, Cylinder* newC,
                                        short oldPosition,
                                        short newPosition,
                                        short motorType,
                                        short boundType,
                                        SubSystem* ps)
    :_oldC(oldC), _newC(newC), _motorType(motorType), _boundType(boundType),
     _oldPosition(oldPosition), _newPosition(newPosition), _ps(ps){}
    
    void operator() (ReactionBase* r) {
        
        //Find motor head
        CCylinder* oldCC = _oldC->getCCylinder();
        CCylinder* newCC = _newC->getCCylinder();
        
        SpeciesMotor* sm1 =oldCC->getCMonomer(_oldPosition)->speciesMotor(_motorType);
        SpeciesMotor* sm2 =newCC->getCMonomer(_newPosition)->speciesMotor(_motorType);
        SpeciesBound* sb2 =newCC->getCMonomer(_newPosition)->speciesBound(_boundType);
        
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        ReactionBase* newOffRxn, *offRxn;
        CCylinder* cc1, *cc2;
        
        //initial set of cylinders
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        
        if(m->getCMotorGhost()->getFirstSpecies() == sm1) {
            
            m->setFirstCylinder(_newC);
            m->setFirstPosition((double)_newPosition
               / SystemParameters::Geometry().cylinderIntSize);
            m->getCMotorGhost()->setFirstSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
            m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
            ({sm2, s1, sb2, s3, s4}, offRxn->getRate());
            
            //remove old off reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        else {
            m->setSecondCylinder(_newC);
            m->setSecondPosition((double)_newPosition
               / SystemParameters::Geometry().cylinderIntSize);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>*)
            m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
            ({s0, sm2, s2, sb2, s4}, offRxn->getRate());
            
            //remove old off reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        
        //set new reaction type
        newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        //attach signal
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
        rcb(newOffRxn->connect(mcallback,false));
        
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        
        //add new
        cc1->addCrossCylinderReaction(cc2, newOffRxn);
        
        //set new unbinding reaction
        m->getCMotorGhost()->setOffReaction(newOffRxn);
        
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
                             SubSystem* ps, Compartment* c = nullptr)
        : _plusEnd(plusEnd), _minusEnd(minusEnd), _filament(filament),
          _compartment(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {

        Compartment* c;
        
        //no compartment was set, pick a random one
        if(_compartment == nullptr)
            c = GController::getRandomCompartment();
        else c = _compartment;
        
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
                SystemParameters::Geometry().cylinderSize, direction);
            
            //check if within boundary
            if(_ps->getBoundary()->within(position) &&
               _ps->getBoundary()->within(npp))
                break;
        }
        
        //create filament, set up ends and filament species
        Filament* f = _ps->addNewFilament(position, direction);
        
        CCylinder* cc = f->getCylinderVector()[0]->getCCylinder();
        int monomerPosition = SystemParameters::Geometry().cylinderIntSize / 2 + 1;
        
        //minus end
        cc->getCMonomer(monomerPosition - 1)->
            speciesMinusEnd(_minusEnd)->up();

        //filament
        cc->getCMonomer(monomerPosition)->
            speciesFilament(_filament)->up();
        cc->getCMonomer(monomerPosition)->
            speciesBound(0)->up();
        
        //plus end
        cc->getCMonomer(monomerPosition + 1)->
            speciesPlusEnd(_plusEnd)->up();
    }
    
};

///Struct to sever a filament based on a reaction
struct FilamentSeveringCallback {
    
    Cylinder* _c1;  ///< Filament severing point
    
    FilamentSeveringCallback(Cylinder* c1) : _c1(c1) {}
    
    void operator() (ReactionBase* r) {
        
        //reactants should be re-marked
        Reaction<SEVERINGREACTANTS + 1, SEVERINGPRODUCTS>* severRxn =
            (Reaction<SEVERINGREACTANTS + 1, SEVERINGPRODUCTS>*)r;
        
        for(int i = 0; i < SEVERINGREACTANTS + 1; i++) severRxn->rspecies()[i]->up();
        
        //sever the filament at given position
        Filament* f = _c1->getFilament();
        Filament* newFilament = f->severFilament(_c1->getPositionFilament());
        
        //if we didn't split, return
        if(newFilament == nullptr) return;
        
        //mark the plus and minus ends of the new and old filament
        CCylinder* cc1 = newFilament->getCylinderVector().back()->getCCylinder();
        CCylinder* cc2 = f->getCylinderVector().front()->getCCylinder();
        
        CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
        CMonomer* m2 = cc2->getCMonomer(0);
        
        short filamentInt1 = m1->activeSpeciesFilament();
        short filamentInt2 = m2->activeSpeciesFilament();
        
        //plus end
        m1->speciesFilament(filamentInt1)->down();
        m1->speciesPlusEnd(filamentInt1)->up();
        m1->speciesBound(0)->down();
        
        //minus end
        m2->speciesFilament(filamentInt2)->down();
        m2->speciesMinusEnd(filamentInt2)->up();
        m2->speciesBound(0)->down();
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
