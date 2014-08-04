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

namespace chem {
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
        ///Initialize the compartment grid, based on the given simulation
        virtual void initializeGrid(CompartmentGrid *grid) = 0;
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CCylinder initialized
        ///@param species - list of species to initialize in CCylinder
        virtual CCylinder* createCCylinder(Compartment* c, CCylinder* lastCCylinder, bool extension) = 0;
        
        ///Remove a CCylinder, based on the given simulation
        virtual void removeCCylinder(CCylinder *cylinder) = 0;

    };
    
    
    ///REACTION CALLBACKS
    
    ///Extension callback
    struct FilamentExtensionCallback {
        
        //members
        Filament* _filament;
        
        ///Constructor, sets members
        FilamentExtensionCallback(Filament* filament) : _filament(filament){};
        
        ///Callback
        void operator() (ReactionBase *r){
            
            //_filament->
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

    
    
} //end namespace chem

#endif /* defined(__Cyto__ChemInitializerImpl__) */
