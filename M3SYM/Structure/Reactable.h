
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

#ifndef M3SYM_Reactable_h
#define M3SYM_Reactable_h

#include "common.h"

/// For a reactable object in the SubSystem.
class Reactable {
    
public:
    ///Update the reactions in this object
    virtual void updateReactionRates() = 0;
};

#endif
