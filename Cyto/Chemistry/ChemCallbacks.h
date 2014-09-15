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
        _filament->printChemComposition();
        
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
        std::cout << "Polymerizing back" << std::endl;
        _filament->PolymerizeBack();
    }
};

///Retraction callback
struct FilamentRetractionCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentRetractionCallback(Filament* filament) : _filament(filament) {};
    
    ///Callback
    void operator() (ReactionBase *r){
        //_filament->
    }
};


///General polymerization callback
struct FilamentPolyCallback {
    
    
    //members
    Filament* _filament;
    
    FilamentPolyCallback(Filament* filament) : _filament(filament) {};
    
    //Callback
    void operator() (ReactionBase *r){
        //_filament->
    }
    
};

///General depolymerization callback
struct FilamentDepolyCallback {
    
    //members
    Filament* _filament;
    
    FilamentDepolyCallback(Filament* filament) : _filament(filament) {};
    
    //Callback
    void operator() (ReactionBase *r){
        //_filament->
    }
    
};


#endif /* defined(__Cyto__Callback__) */
