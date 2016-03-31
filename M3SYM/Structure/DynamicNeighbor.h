//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_DynamicNeighbor_h
#define M3SYM_DynamicNeighbor_h

#include "common.h"

#include "Neighbor.h"

/// An abstract base class for any element that can be
/// added or removed from a NeighborList dynamically at runtime.
class DynamicNeighbor : public Neighbor {
    
protected:
    DynamicNeighbor() : Neighbor() {}

public:
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~DynamicNeighbor() noexcept {}
};

#endif
