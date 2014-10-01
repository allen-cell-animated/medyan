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

///REACTION CALLBACKS

///Extension callback
struct FilamentExtensionFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentExtensionFrontCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        _filament->PolymerizeFront();
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
        _filament->PolymerizeBack();
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
        _filament->DepolymerizeFront();
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
        _filament->DepolymerizeBack();
    }
};



#endif /* defined(__Cyto__Callback__) */
