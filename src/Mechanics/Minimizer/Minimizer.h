
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Minimizer_h
#define MEDYAN_Minimizer_h

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
    /// @param stepLimit - If there is a limit for the number of equlibration steps.
    void virtual equlibrate(ForceFieldManager &FFM, bool stepLimit) = 0;

	tuple<floatingpoint, vector<floatingpoint>, vector<string>>  virtual getEnergy(ForceFieldManager &FFM, floatingpoint d) = 0;
    
    
    
};

#endif
