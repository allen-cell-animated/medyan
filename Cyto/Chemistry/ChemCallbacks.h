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
#include "common.h"
#include "Filament.h"
#include "ReactionBase.h"
#include "SubSystem.h"

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
    Cylinder* _c1;
    Cylinder* _c2;
    short _linkerType;

    LinkerBindingCallback(SubSystem* ps, Cylinder* c1, Cylinder* c2, short linkerType)
        : _ps(ps), _c1(c1), _c2(c2), _linkerType(linkerType) {}
    
    void operator() (ReactionBase *r) {
        
        //find locations of two species
        
        
        
        
        
        
        
    }
    
    
};










#endif /* defined(__Cyto__Callback__) */
