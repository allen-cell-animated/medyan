
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

#ifndef M3SYM_Movable_h
#define M3SYM_Movable_h

#include "common.h"

/// For a movable object in the SubSystem.
class Movable {
    
public:
    /// Update the position of this object
    virtual void updatePosition() = 0;
};

#endif
