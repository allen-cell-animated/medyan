
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

#ifndef M3SYM_Trackable_h
#define M3SYM_Trackable_h

#include "common.h"

/// An abstract base class for a trackable object in the SubSystem.

/*!
 *  Every class extending Trackable will provide a container to track, 
 *  add, and remove instances of its class from the SubSystem. In general,
 *  this is done using the Database class.
 */

class Trackable {
    
protected:
    Trackable() {}
    
public:
    //@{
    /// Add or remove this Trackable element from the SubSystem.
    /// @note - this update could be due to an updated force minimization,
    /// a set of chemical steps, or any other event in the SubSystem.
    virtual void addToSubSystem() = 0;
    
    virtual void removeFromSubSystem() = 0;
    //@}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Trackable() noexcept {}
};


#endif
