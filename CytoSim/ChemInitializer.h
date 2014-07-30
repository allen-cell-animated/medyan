//
//  ChemInitializer.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemInitializer__
#define __Cyto__ChemInitializer__

#include <iostream>
#include <vector>

namespace chem {
    
    ///Forward declarations
    class ChemInitializerImpl;

    class Compartment;
    class CCylinder;
    
    ///Key for initialization of ChemInitializer
    class ChemInitializerInitKey { friend class ChemController;
                                    ChemInitializerInitKey(){}
                                    ~ChemInitializerInitKey(){} };
    
    ///Key for the initialization of grid
    class ChemInitializerGridKey { friend class ChemController;
                                    ChemInitializerGridKey(){}
                                    ~ChemInitializerGridKey(){} };
    
    ///Key for the creation and destruction of CCylinders
    class ChemInitializerCylinderKey { friend class Cylinder;
                                        ChemInitializerCylinderKey(){}
                                        ~ChemInitializerCylinderKey(){} };
    
    /// ChemInitializer class is used for initailizing chemical reactions based on a specific system
    /*!
     *  CFilamentInitializer is a singleton used for initailizing all chemistry in 
     *  the system, including CCylinders as well as the compartment grid. Like ChemSim, specific functions 
     *  in this class require a key to be accessed. The key declarations are above. These keys can
     *  only be created and/or destroyed by the classes that are "friends" with the key.
     */
    class ChemInitializer {
        
    public:
        ///Set the chemInitializer instance
        void setInstance(ChemInitializerInitKey k, ChemInitializerImpl *cii);
        
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CCylinderinitialized
        ///@param species - list of species to initialize in CCylinder
        CCylinder* createCCylinder(ChemInitializerCylinderKey k,
                                   Compartment* c,
                                   std::vector<std::string> species,
                                   int length);
        
        ///Remove a CCylinder, based on the given simulation
        void removeCCylinder(ChemInitializerCylinderKey k, CCylinder *cylinder);
        
        
    private:
        static ChemInitializerImpl* _pimpl; ///< Store a pointer to a specific implementation of the initializer; no ownership
        ChemInitializer() {};

    };
    
    ///REACTION CALLBACKS
    
    ///Extension callback
    template<size_t NDIM>
    struct FilamentExtensionCallback {
        
        //members
        Filament* _filament;
        std::vector<std::string> _species;
        
        ///Constructor, sets members
        CFilamentExtensionCallback(Filament* filament,
                                   std::vector<std::string> species) :
       _filament(filament), _species(species) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            
            _filament->
        }
    };
    
    ///Retraction callback
    template<size_t NDIM>
    struct FilamentRetractionCallback {
        
        //members
        Filament* _filament;
        
        ///Constructor, sets members
        CFilamentRetractionCallback(Filament* filament) :
        _filament(filament) {};
        
        ///Callback
        void operator() (ReactionBase *r){
            _filament->
        }
    };
    
    
    ///General polymerization callback
    template<size_t NDIM>
    struct FilamentPolyCallback {
        
        //members
        CFilament* _filament;
        
        CFilamentPolyCallback(CFilament* filament) :
        _filament(filament) {};
        
        //Callback
        void operator() (ReactionBase *r){
           _filament->
        }
        
    };
    
    ///General depolymerization callback
    template<size_t NDIM>
    struct FilamentDepolyCallback {
        
        //members
        CFilament* _filament;
        
        CFilamentDepolyCallback(CSystem<NDIM>* csystem,
                                CFilament* filament) :
        _filament(filament) {};
        
        //Callback
        void operator() (ReactionBase *r){
            _filament->
        }
        
    };

    
    
    
    
}; //end namespace chem

#endif /* defined(__Cyto__ChemInitializer__) */