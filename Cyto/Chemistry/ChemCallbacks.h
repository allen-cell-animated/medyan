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
        std::cout << "Polymerizing front" << std::endl;
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
        _filament->printChemComposition();
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
        std::cout << "Depolymerizing front" << std::endl;
        _filament->DepolymerizeFront();
        _filament->printChemComposition();
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
        std::cout << "Depolymerizing back" << std::endl;
        _filament->DepolymerizeBack();
        _filament->printChemComposition();
    }
};


#endif /* defined(__Cyto__Callback__) */
