
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

/// Callback to bind a Linker to Filament
struct LinkerBindingCallback {
    
    SubSystem* _ps;
    CCylinder* _cc1, *_cc2;
    short _linkerType;
    short _position1, _position2;

    LinkerBindingCallback(CCylinder* cc1, CCylinder* cc2, short linkerType, short position1, short position2, SubSystem* ps)
        : _ps(ps), _cc1(cc1), _cc2(cc2), _linkerType(linkerType),  _position1(position1), _position2(position2){}
    
    void operator() (ReactionBase *r) {
        
        // Create a linker
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        _ps->addNewLinker(_cc1->getCylinder(), _cc2->getCylinder(), _linkerType, pos1, pos2);
    }
};

/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;
    CCylinder* _cc1, *_cc2;
    short _motorType;
    short _position1, _position2;
    
    MotorBindingCallback(CCylinder* cc1, CCylinder* cc2, short motorType, short position1, short position2, SubSystem* ps)
    : _ps(ps), _cc1(cc1), _cc2(cc2), _motorType(motorType),  _position1(position1), _position2(position2){}
    
    void operator() (ReactionBase *r) {
        
        // Create a motor
        int cylinderSize = SystemParameters::Geometry().cylinderIntSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        _ps->addNewMotorGhost(_cc1->getCylinder(), _cc2->getCylinder(), _motorType, pos1, pos2);
    }
};


/// Callback to unbind a CBound from a Filament
struct UnbindingCallback {
    
    SubSystem* _ps;
    SpeciesBound* _s1;
    
    UnbindingCallback(SpeciesBound* s1, SubSystem* ps) : _s1(s1), _ps(ps) {}
    
    void operator() (ReactionBase *r) {
        
        //check if we have a basic bound element, linker, or motor
        CBound* cBound = _s1->getCBound();
        
        if(cBound != nullptr) {
            if(dynamic_cast<CLinker*>(cBound))
                _ps->removeLinker(static_cast<CLinker*>(cBound)->getLinker());

            else if(dynamic_cast<CMotorGhost*>(cBound))
                _ps->removeMotorGhost(static_cast<CMotorGhost*>(cBound)->getMotorGhost());
        }
    }
};


/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingForwardCallback {
    
    SpeciesMotor* _sm1;
    SpeciesMotor* _sm2;
    
    MotorWalkingForwardCallback(SpeciesMotor* sm1, SpeciesMotor* sm2)
        :_sm1(sm1), _sm2(sm2) {}
    
    void operator() (ReactionBase* r) {
        
        if(_sm1->getCBound() == nullptr) {
            
            cout << "Major bug: motor is in wrong place" << endl;
            return;
        }
        
        MotorGhost* m = static_cast<CMotorGhost*>(_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double shift = 1.0 / SystemParameters::Geometry().cylinderIntSize;
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

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingBackwardCallback {
    
    SpeciesMotor* _sm1;
    SpeciesMotor* _sm2;
    
    MotorWalkingBackwardCallback(SpeciesMotor* sm1, SpeciesMotor* sm2)
        :_sm1(sm1), _sm2(sm2) {}
    
    void operator() (ReactionBase* r) {
        
        if(_sm1->getCBound() == nullptr) {
            
            cout << "Major bug: motor is in wrong place" << endl;
            return;
        }
        
        MotorGhost* m = static_cast<CMotorGhost*>(_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double shift = 1.0 / SystemParameters::Geometry().cylinderIntSize;
        double newPosition;
        
        if(m->getCMotorGhost()->getFirstSpecies() == _sm1) {
            newPosition = m->getFirstPosition() - shift;
            m->setFirstPosition(newPosition);
            m->getCMotorGhost()->setFirstSpecies(_sm2);
        }
        else {
            newPosition = m->getSecondPosition() - shift;
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
    
    MotorMovingCylinderForwardCallback(SpeciesMotor* sm1, SpeciesMotor* sm2, CCylinder* newCCylinder)
        : _sm1(sm1), _sm2(sm2), _newCCylinder(newCCylinder) {}
    
    void operator() (ReactionBase* r) {
        
        if(_sm1->getCBound() == nullptr) {
            
            cout << "Major bug: motor is in wrong place" << endl;
            return;
        }
        
        MotorGhost* m = static_cast<CMotorGhost*>(_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double newPosition = 0.0;
        
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


/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderBackwardCallback {
    
    //members
    SpeciesMotor* _sm1;
    SpeciesMotor* _sm2;
    CCylinder* _newCCylinder;
    
    MotorMovingCylinderBackwardCallback(SpeciesMotor* sm1, SpeciesMotor* sm2, CCylinder* newCCylinder)
        : _sm1(sm1), _sm2(sm2), _newCCylinder(newCCylinder) {}
    
    void operator() (ReactionBase* r) {
        
        if(_sm1->getCBound() == nullptr) {
            
            cout << "Major bug: motor is in wrong place" << endl;
            return;
        }
        
        MotorGhost* m = static_cast<CMotorGhost*>(_sm1->getCBound())->getMotorGhost();
        
        //shift the position of one side of the motor forward
        double newPosition = 1.0 - 1.0 / SystemParameters::Geometry().cylinderIntSize;
        
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
