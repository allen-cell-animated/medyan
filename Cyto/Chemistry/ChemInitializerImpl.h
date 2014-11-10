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
struct ChemistryData;

///ChemInitializerImpl is an abstract base class for initialization of all chemistry in the system

/*  
 *  Specific initializers should inherit from ChemInitializerImpl. A user will then attach the corresponding 
 *  initializer to ChemInitializer via the initializer base class, ChemInitializerImpl.
 */
class ChemInitializerImpl {
    
public:
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~ChemInitializerImpl() noexcept{}
    
    ///Initialize the compartment grid, based on the given simulation
    virtual void initialize(ChemistryData& chem) = 0;
    
    ///Initializer, based on the given simulation
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c,
                                       bool extensionFront, bool extensionBack, bool creation) = 0;

    ///add/update cross cylinder reactions that are within range
    virtual void updateCCylinder(CCylinder* cc, vector<CCylinder*>& cNeighbors) = 0;
    
    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};

#endif /* defined(__Cyto__ChemInitializerImpl__) */
