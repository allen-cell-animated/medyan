
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
    
    vector<Species*> _offSpecies;  ///< Species to add to the unbinding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    LinkerBindingCallback(Cylinder* c1, Cylinder* c2, short linkerType, short position1, short position2,
                          vector<Species*> offSpecies, float offRate, SubSystem* ps)
                          : _ps(ps), _c1(c1), _c2(c2), _linkerType(linkerType),
                            _position1(position1), _position2(position2),
                            _offSpecies(offSpecies), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
        
        // Create a linker
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        Linker* l = _ps->addNewLinker(_c1, _c2, _linkerType, pos1, pos2);
        
        //Attach the callback to the off reaction, add it
        ReactionBase* offRxn = new Reaction<2,3>(_offSpecies, _offRate);
        offRxn->setReactionType(ReactionType::LINKERUNBINDING);
        
        LinkerUnbindingCallback lcallback(l, _ps);
        boost::signals2::shared_connection_block rcb(offRxn->connect(lcallback,false));
        
        _c1->getCCylinder()->addCrossCylinderReaction(_c2->getCCylinder(), offRxn);
    }
};

/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;
    CCylinder* _cc1, *_cc2;
    short _motorType;
    short _position1, _position2;
    
    MotorBindingCallback(CCylinder* cc1, CCylinder* cc2, short motorType,
                         short position1, short position2, SubSystem* ps)
                         : _ps(ps), _cc1(cc1), _cc2(cc2), _motorType(motorType),
                           _position1(position1), _position2(position2){}
    
    void operator() (ReactionBase *r) {
        
        // Create a motor
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        _ps->addNewMotorGhost(_cc1->getCylinder(), _cc2->getCylinder(), _motorType, pos1, pos2);
    }
};


/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    SpeciesMotor* _s1;
    
    MotorUnbindingCallback(SpeciesMotor* s1, SubSystem* ps) : _s1(s1), _ps(ps) {}
    
    void operator() (ReactionBase *r) {
        _ps->removeMotorGhost(((CMotorGhost*)_s1->getCBound())->getMotorGhost());
    }
};

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingForwardCallback {
    
    SpeciesMotor* _sm1;
    SpeciesMotor* _sm2;
    
    MotorWalkingForwardCallback(SpeciesMotor* sm1, SpeciesMotor* sm2)
                                :_sm1(sm1), _sm2(sm2) {}
    
    void operator() (ReactionBase* r) {
        
        MotorGhost* m = ((CMotorGhost*)_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double shift =  1.0 / SystemParameters::Chemistry().numBindingSites;
        double newPosition;
        
        if(m->getCMotorGhost()->getFirstSpecies() == _sm1) {
            newPosition = m->getFirstPosition() + shift;
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(_sm2);
        }
        else {
            newPosition = m->getSecondPosition() + shift;
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(_sm2);
        }
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderForwardCallback {
    
    //members
    SpeciesMotor* _sm1;
    SpeciesMotor* _sm2;
    CCylinder* _newCCylinder;
    
    MotorMovingCylinderForwardCallback(SpeciesMotor* sm1,
                                       SpeciesMotor* sm2,
                                       CCylinder* newCCylinder)
                                        : _sm1(sm1), _sm2(sm2), _newCCylinder(newCCylinder) {}
    
    void operator() (ReactionBase* r) {
        
        MotorGhost* m = ((CMotorGhost*)_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double newPosition = 1.0 / (2 * SystemParameters::Chemistry().numBindingSites);
        
        if(m->getCMotorGhost()->getFirstSpecies() == _sm1) {
            m->setFirstCylinder(_newCCylinder->getCylinder());
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(_sm2);
        }
        else {
            m->setSecondCylinder(_newCCylinder->getCylinder());
            m->setSecondPosition(newPosition);
            m->getCMotorGhost()->setSecondSpecies(_sm2);
        }
    }
};

#endif
