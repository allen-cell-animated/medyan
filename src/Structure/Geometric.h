
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Geometric_h
#define MEDYAN_Geometric_h

#include "common.h"

/// An abstract base class for a geometric element in the SubSystem.

/*! The main function of the Geometric class is to implement updateGeometry(),
 *  so that the geometry of any object extending this class can be updated
 *  by the SubSystem.
 * 
 *  Note that the geometry of some objects may depend on others, so it is
 *  preferable to inherit Geometric only for some parental objects that know
 *  clearly the order of geometry updating for their childern.
 */
class Geometric {
    
protected:
    Geometric() {}
    
public:
    /// Update the geometry of this element
    /// @note - many processes could change the geometric properties of the
    /// object, but its main purpose is to provide intermediate calculations
    /// for mechanical relaxation.
    /// Function will be called by the SubSystem on all Geometric objects.
    virtual void updateGeometry(bool calcDerivative=false, double d=0.0) = 0;
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Geometric() noexcept {}
};

#endif
