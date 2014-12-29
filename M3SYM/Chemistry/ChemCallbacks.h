
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

#include "SubSystem.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"

#include "SystemParameters.h"

/// Callback to extend the front of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionFrontCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentExtensionFrontCallback(Filament* filament) : _filament(filament){};
    //Callback
    void operator() (ReactionBase *r){ _filament->extendFront(); }
};

/// Callback to extend the back of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionBackCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentExtensionBackCallback(Filament* filament) : _filament(filament){};
    //Callback
    void operator() (ReactionBase *r){ _filament->extendBack(); }
};

/// Callback to retract the front of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionFrontCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentRetractionFrontCallback(Filament* filament) : _filament(filament) {};
    //Callback
    void operator() (ReactionBase *r){ _filament->retractFront(); }
};

/// Callback to retract the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionBackCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentRetractionBackCallback(Filament* filament) : _filament(filament) {};
    //Callback
    void operator() (ReactionBase *r){ _filament->retractBack(); }
};

/// Callback to polymerize the front of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationFrontCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentPolymerizationFrontCallback(Filament* filament) : _filament(filament){};
    //Callback
    void operator() (ReactionBase *r){ _filament->polymerizeFront();}
};

/// Callback to polymerize the back of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationBackCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentPolymerizationBackCallback(Filament* filament) : _filament(filament){};
    //Callback
    void operator() (ReactionBase *r){ _filament->polymerizeBack(); }
};

/// Callback to depolymerize the front of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationFrontCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentDepolymerizationFrontCallback(Filament* filament) : _filament(filament) {};
    //Callback
    void operator() (ReactionBase *r){ _filament->depolymerizeFront(); }
};

/// Callback to depolymerize the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationBackCallback {
    
    Filament* _filament;
    
    //Constructor, sets members
    FilamentDepolymerizationBackCallback(Filament* filament) : _filament(filament) {};
    //Callback
    void operator() (ReactionBase *r){ _filament->depolymerizeBack(); }
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
        
        //reset the associated reactions
        _linker->getMLinker()->stretchForce = 0.0;
        _linker->updateReactionRates();
        
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
        : _ps(ps), _c1(c1), _c2(c2),
          _linkerType(linkerType),
          _position1(position1), _position2(position2), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        // Create a linker
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        Linker* l = _ps->addNewLinker(_c1, _c2, _linkerType, pos1, pos2);
        
        //Attach the callback to the off reaction, add it
        //first, get species from this rxn
        Reaction<3,2>* onRxn = (Reaction<3,2>*)r;
        RSpecies** rSpecies = onRxn->rspecies();
        vector<Species*> offSpecies;
        
        //copy into offspecies vector
        for(int i = 3; i < 5; i++) offSpecies.push_back(&rSpecies[i]->getSpecies());
        for(int i = 0; i < 3; i++) offSpecies.push_back(&rSpecies[i]->getSpecies());
        
        ReactionBase* offRxn = new Reaction<2,3>(offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::LINKERUNBINDING);
        
        LinkerUnbindingCallback lcallback(l, _ps);
        boost::signals2::shared_connection_block rcb(offRxn->connect(lcallback,false));
        
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
        
        //reset the associated reactions
        _motor->getMMotorGhost()->stretchForce = 0.0;
        _motor->updateReactionRates();
        
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
        Reaction<3,2>* onRxn = (Reaction<3,2>*)r;
        RSpecies** rSpecies = onRxn->rspecies();
        vector<Species*> offSpecies;
        
        //copy into offspecies vector
        for(int i = 3; i < 5; i++) offSpecies.push_back(&(rSpecies[i]->getSpecies()));
        for(int i = 0; i < 3; i++) offSpecies.push_back(&(rSpecies[i]->getSpecies()));
        
        ReactionBase* offRxn = new Reaction<2,3>(offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::MOTORUNBINDING);
        
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block rcb(offRxn->connect(mcallback,false));
        
        _c1->getCCylinder()->addCrossCylinderReaction(_c2->getCCylinder(), offRxn);
        m->getCMotorGhost()->setOffReaction(offRxn);
        
        //update the reaction rates
        m->updateReactionRates();
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
        SpeciesMotor* sm1 =
        _c->getCCylinder()->getCMonomer(_oldPosition)->speciesMotor(_motorType);
        SpeciesMotor* sm2 =
        _c->getCCylinder()->getCMonomer(_newPosition)->speciesMotor(_motorType);
        SpeciesBound* sb2 =
        _c->getCCylinder()->getCMonomer(_newPosition)->speciesBound(_boundType);
        
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
            offRxn = (Reaction<2,3>*)m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn =
            new Reaction<2,3>({sm2, s1, sb2, s3, s4}, offRxn->getBareRate());
        }
        else {
            newPosition = m->getSecondPosition() + shift;
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<2,3>*)m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn =
            new Reaction<2,3>({s0, sm2, s2, sb2, s4}, offRxn->getBareRate());
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
        
        //update reaction rates
        m->updateReactionRates();
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderForwardCallback {
    
    Cylinder* _oldC;        ///< Old cylinder the motor is attached to
    Cylinder* _newC;        ///< New cylinder motor will be attached to
    
    short _oldPosition;     ///< Old position of motor head
    
    short _motorType;       ///< Type of motor
    short _boundType;       ///< Type of bound this motor is taking place of
    
    SubSystem* _ps;         ///< Ptr to subsystem
    
    MotorMovingCylinderForwardCallback(Cylinder* oldC, Cylinder* newC,
                                       short oldPosition,
                                       short motorType, short boundType, SubSystem* ps)
        :_oldC(oldC), _newC(newC), _motorType(motorType),
         _boundType(boundType), _oldPosition(oldPosition), _ps(ps){}
    
    void operator() (ReactionBase* r) {
        
        //Find the motor
        double newPosition = 1.0 / (2 * SystemParameters::Chemistry().numBindingSites);
        int newIntPosition =
            newPosition * SystemParameters::Geometry().cylinderIntSize + 1;

        SpeciesMotor* sm1 =
            _oldC->getCCylinder()->getCMonomer(_oldPosition)->
            speciesMotor(_motorType);
        SpeciesMotor* sm2 =
            _newC->getCCylinder()->getCMonomer(newIntPosition)->
            speciesMotor(_motorType);
        SpeciesBound* sb2 =
            _newC->getCCylinder()->getCMonomer(newIntPosition)->
            speciesBound(_boundType);
        
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        ReactionBase* newOffRxn, *offRxn;
        CCylinder* cc1, *cc2;
        //initial set of cylinders
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        
        if(m->getCMotorGhost()->getFirstSpecies() == sm1) {
            
            m->setFirstCylinder(_newC);
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<2,3>*)m->getCMotorGhost()->getOffReaction();
            
            //take the first, third, and fourth species
            Species* s1 = &(offRxn->rspecies()[1]->getSpecies());
            Species* s3 = &(offRxn->rspecies()[3]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn =
                new Reaction<2,3>({sm2, s1, sb2, s3, s4}, offRxn->getBareRate());
            
            //remove old reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        else {
            m->setSecondCylinder(_newC);
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(sm2);
            
            //change off reaction to include new species
            offRxn = (Reaction<2,3>*)m->getCMotorGhost()->getOffReaction();
            
            //take the zeroth, second, and fourth species
            Species* s0 = &(offRxn->rspecies()[0]->getSpecies());
            Species* s2 = &(offRxn->rspecies()[2]->getSpecies());
            Species* s4 = &(offRxn->rspecies()[4]->getSpecies());
            
            //create new reaction
            newOffRxn =
                new Reaction<2,3>({s0, sm2, s2, sb2, s4}, offRxn->getBareRate());
            
            //remove old reaction
            cc1->removeCrossCylinderReaction(cc2, offRxn);
        }
        //set new reaction type
        newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
        //attach signal
        MotorUnbindingCallback mcallback(m, _ps);
        boost::signals2::shared_connection_block
            rcb(newOffRxn->connect(mcallback,false));
        
        //add new reaction
        cc1 = m->getFirstCylinder()->getCCylinder();
        cc2 = m->getSecondCylinder()->getCCylinder();
        cc1->addCrossCylinderReaction(cc2, newOffRxn);
        
        //set new unbinding reaction
        m->getCMotorGhost()->setOffReaction(newOffRxn);
        
        //update reaction rates
        m->updateReactionRates();
    }
};

#endif
