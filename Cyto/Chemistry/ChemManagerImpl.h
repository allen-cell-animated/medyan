//
//  ChemManagerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemManagerImpl__
#define __Cyto__ChemManagerImpl__

#include <iostream>
#include <vector>

#include "common.h"
#include "CompartmentContainer.h"

class Filament;
class CCylinder;
class ReactionBase;
struct ChemistryData;

///ChemManagerImpl is an abstract base class for initialization of all chemistry in the system

/*  
 *  Specific Managers should inherit from ChemManagerImpl. A user will then attach the corresponding 
 *  Manager to ChemManager via the Manager base class, ChemManagerImpl.
 */
class ChemManagerImpl {
    
public:
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~ChemManagerImpl() noexcept{}
    
    ///Initialize the compartment grid, based on the given simulation
    virtual void initialize(ChemistryData& chem) = 0;
    
    ///Manager, based on the given simulation
    virtual void  initializeCCylinder(CCylinder* cc, Filament* f,
                                      bool extensionFront, bool extensionBack, bool creation) = 0;

    ///add/update cross cylinder reactions that are within range
    virtual void updateCCylinder(CCylinder* cc) = 0;
    
    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};

#endif /* defined(__Cyto__ChemManagerImpl__) */
