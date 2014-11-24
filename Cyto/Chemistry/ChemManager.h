//
//  ChemManager.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemManager__
#define __Cyto__ChemManager__

#include <iostream>

#include "common.h"

#include "ChemManagerImpl.h"

///FORWARD DECLARATIONS
class Compartment;
class Filament;
class CCylinder;
struct ChemistryData;

/// ChemManager class is used for initailizing chemical reactions based on a specific system
/*!
 *  ChemManager is a singleton used for initailizing all chemistry in
 *  the system, including CCylinders as well as the compartment grid. Like ChemSim, specific functions 
 *  in this class require a key to be accessed. The key declarations are above. These keys can
 *  only be created and/or destroyed by the classes that are "friends" with the key.
 */
class ChemManager {
    
public:
    ///Set the chemManager instance
    static void setInstance(ChemManagerImpl *cii);
    
    ///Initialize the compartment grid, based on the given simulation
    static void initialize(ChemistryData& chem);
    
    ///Initializer, based on the given simulation
    static void initializeCCylinder(CCylinder* cc, Filament* f,bool extensionFront, bool extensionBack, bool creation);
    
    ///add/update cross cylinder reactions that are within range
    static void updateCCylinder(CCylinder* cc);
    
    
private:
    static ChemManagerImpl* _pimpl; ///< Store a pointer to a specific implementation of the initializer; no ownership
    ChemManager() {};

};


#endif /* defined(__Cyto__ChemManager__) */