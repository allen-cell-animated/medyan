//
//  ChemCallbacks.h
//  Cyto
//
//  Created by James Komianos on 9/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemCallbacks__
#define __Cyto__ChemCallbacks__

#include <iostream>
#include "ReactionBase.h"
#include "SubSystem.h"

#include "common.h"
#include "SystemParameters.h"
using namespace std;

///FILAMENT REACTION CALLBACKS

///Extension callback
struct FilamentExtensionFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentExtensionFrontCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->ExtendFront();
    }
};

///Extension callback
struct FilamentExtensionBackCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentExtensionBackCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->ExtendBack();
    }
};

///Retraction callback
struct FilamentRetractionFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentRetractionFrontCallback(Filament* filament) : _filament(filament) {};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->RetractFront();
    }
};

///Retraction callback
struct FilamentRetractionBackCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentRetractionBackCallback(Filament* filament) : _filament(filament) {};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->RetractBack();
    }
};

///Polymerization/depolymerization callbacks
struct FilamentPolymerizationFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentPolymerizationFrontCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->PolymerizeFront();
    }
};

struct FilamentPolymerizationBackCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentPolymerizationBackCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->PolymerizeBack();
    }
};

///Retraction callback
struct FilamentDepolymerizationFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentDepolymerizationFrontCallback(Filament* filament) : _filament(filament) {};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->DepolymerizeFront();
    }
};

///Retraction callback
struct FilamentDepolymerizationBackCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentDepolymerizationBackCallback(Filament* filament) : _filament(filament) {};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->DepolymerizeBack();
    }
};

///LINKER AND MOTOR CALLBACKS

///Linker binding callback
struct LinkerBindingCallback {
    
    ///members
    SubSystem* _ps;
    Cylinder* _c1, *_c2;
    short _linkerType;
    short _position1, _position2;

    LinkerBindingCallback(Cylinder* c1, Cylinder* c2, short linkerType, short position1, short position2, SubSystem* ps)
        : _ps(ps), _c1(c1), _c2(c2), _linkerType(linkerType),  _position1(position1), _position2(position2){}
    
    void operator() (ReactionBase *r) {
        
        ///Create a linker
        int cylinderSize = SystemParameters::Geometry().cylinderSize /
                           SystemParameters::Geometry().monomerSize;
        
        double pos1 = double(_position1) / cylinderSize;
        double pos2 = double(_position2) / cylinderSize;
        
        _ps->AddNewLinker(_c1, _c2, _linkerType, pos1, pos2);
        
        Linker* newLinker = LinkerDB::Instance(LinkerDBKey())->back();
        
        ///attach species
        //SpeciesBound* s1 = r->
        
        
        
        
    }
    
    
};










#endif /* defined(__Cyto__Callback__) */
