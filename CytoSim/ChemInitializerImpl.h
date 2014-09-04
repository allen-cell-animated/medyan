//
//  ChemInitializerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemInitializerImpl__
#define __Cyto__ChemInitializerImpl__

#include <iostream>
#include <vector>
#include "MFilament.h"
#include "CompartmentContainer.h"

class ReactionBase;
class CompartmentGrid;
class Compartment;
class CCylinder;

///ChemInitializerImpl is an abstract base class for initialization of all chemistry in the system

/*  
 *  Specific initializers should inherit from ChemInitializerImpl. A user will then attach the corresponding 
 *  initializer to ChemInitializer via the initializer base class, ChemInitializerImpl.
 */
class ChemInitializerImpl {
    
public:
    virtual ~ChemInitializerImpl() {}
    
    ///Initialize the compartment grid, based on the given simulation
    virtual void initializeGrid() = 0;
    
    ///Initializer, based on the given simulation
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack) = 0;
    
    ///Remove a CCylinder, based on the given simulation
    virtual void removeCCylinder(Filament* pf, bool retractionFront, bool retractionBack) = 0;

    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};


///REACTION CALLBACKS

///Extension callback
struct FilamentExtensionFrontCallback {
    
    //members
    Filament* _filament;
    
    ///Constructor, sets members
    FilamentExtensionFrontCallback(Filament* filament) : _filament(filament){};
    
    ///Callback
    void operator() (ReactionBase *r){
        //_filament->PolymerizeFront();
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
        //_filament->PolymerizeBack();
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

#endif /* defined(__Cyto__ChemInitializerImpl__) */
