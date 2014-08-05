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

///Forward declarations
class ChemInitializerImpl;
class CompartmentGrid;
class Compartment;
class CCylinder;

///Key for initialization of ChemInitializer
class ChemInitializerInitKey { friend class CController; ChemInitializerInitKey(){}; ~ChemInitializerInitKey(){}; };

///Key for the initialization of grid
class ChemInitializerGridKey { friend class CController; ChemInitializerGridKey(){}; ~ChemInitializerGridKey(){}; };

///Key for the creation and destruction of CCylinders
class ChemInitializerCylinderKey { friend class Filament; ChemInitializerCylinderKey(){}; ~ChemInitializerCylinderKey(){}; };

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
    static void setInstance(ChemInitializerInitKey k, ChemInitializerImpl *cii);
    
    ///Initialize the compartment grid, based on the given simulation
    static void initializeGrid(ChemInitializerGridKey k);
    
    ///Initializer, based on the given simulation
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    static CCylinder* createCCylinder(ChemInitializerCylinderKey k,
                               Compartment* c,
                               CCylinder* lastCCylinder, bool extension = false);
    
    ///Remove a CCylinder, based on the given simulation
    static void removeCCylinder(ChemInitializerCylinderKey k, CCylinder *cylinder);
    
    
private:
    static ChemInitializerImpl* _pimpl; ///< Store a pointer to a specific implementation of the initializer; no ownership
    ChemInitializer() {};

};


#endif /* defined(__Cyto__ChemInitializer__) */