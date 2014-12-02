
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Neighbor_h
#define M3SYM_Neighbor_h

#include "common.h"

/// An abstract class for any object that can be added/removed from a [NeighborList](@ref NeighborList).
/// @note - any subclass that inherits from Neighbor MUST add and remove itself to the neighbors list DB.
class Neighbor {
    
protected:
    Neighbor() {}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~Neighbor() noexcept {}
};

#endif
