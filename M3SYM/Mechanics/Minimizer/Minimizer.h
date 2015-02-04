
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

#ifndef M3SYM_Minimizer_h
#define M3SYM_Minimizer_h

#include "common.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

/// A mechanical minimzer used by the [MController](@ref MController).
/*!
 *  This class is used to mechanically equilibrate a system. Implementations of this 
 *  class can be either conjugate gradient methods, or Langevin Dynamics minimizers.
 */
class Minimizer {
public:
    /// Equilibrate a system based on given configuration and force-fields
    /// @param FFM - The force field manager of this system. @see ForceFieldManager.
    void virtual equlibrate(ForceFieldManager &FFM) = 0;
};

#endif
