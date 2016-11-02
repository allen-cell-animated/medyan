
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Reactable_h
#define MEDYAN_Reactable_h

#include "common.h"

/// An abstract base class for a reactable element in the SubSystem.

/*! The main function of the Reactable class is to implement updateReactionRates(),
 *  so that the reactions related to any object extending this class can be updated
 *  by the SubSystem.
 */
class Reactable {
    
protected:
    Reactable() {}
    
public:
    ///Update the reactions in this element
    /// @note - this update could be due to an updated force minimization,
    /// a set of chemical steps, or any other event in the SubSystem. This
    /// function will be called by the SubSystem on all Reactables.
    virtual void updateReactionRates() = 0;
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Reactable() noexcept {}
};

#endif
