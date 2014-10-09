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

#include "common.h"
#include "CompartmentContainer.h"

class Filament;
class CCylinder;
class ReactionBase;
struct ChemistrySpeciesAndReactions;

///ChemInitializerImpl is an abstract base class for initialization of all chemistry in the system

/*  
 *  Specific initializers should inherit from ChemInitializerImpl. A user will then attach the corresponding 
 *  initializer to ChemInitializer via the initializer base class, ChemInitializerImpl.
 */
class ChemInitializerImpl {
    
public:
    virtual ~ChemInitializerImpl() {}
    
    ///Initialize the compartment grid, based on the given simulation
    virtual void initialize(ChemistrySpeciesAndReactions& chemSR) = 0;
    
    ///Initializer, based on the given simulation
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack) = 0;

    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};

#endif /* defined(__Cyto__ChemInitializerImpl__) */
