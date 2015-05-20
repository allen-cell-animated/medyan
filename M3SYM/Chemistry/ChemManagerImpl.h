
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_ChemManagerImpl_h
#define M3SYM_ChemManagerImpl_h

#include "common.h"

#include "CompartmentContainer.h"

//FORWARD DECLARATIONS
class Filament;
class CCylinder;
class ReactionBase;
struct ChemistryData;

/// An abstract base class for initialization of all chemistry in the system
/*  
 *  Specific Managers should inherit from ChemManagerImpl. A user will then attach the
 *  corresponding Manager to ChemManager via the Manager base class, ChemManagerImpl.
 *  @note @see ChemManager for documentation on the chem manager singleton.
 */
class ChemManagerImpl {
    
public:
    /// Destructor does nothing.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ChemManagerImpl() noexcept{}
    
    virtual void initializeSystem() = 0;
    
    virtual void initializeCCylinder(CCylinder* cc, Filament* f,
                                     bool extensionFront,
                                     bool extensionBack,
                                     bool initialization) = 0;
    
    virtual void updateCopyNumbers() = 0;
};

#endif 
