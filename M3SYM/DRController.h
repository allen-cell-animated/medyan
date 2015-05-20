
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

#ifndef M3SYM_DRController_h
#define M3SYM_DRController_h

#include "common.h"

//FORWARD DECLARATIONS
struct DynamicRateTypes;

/// Used to initialize the dynamic rate components of a simulation.

/*!
 *  The DRController intializes all dynamic rate objects in the system, 
 *  as specified in the input file. The controller's intiailize function will 
 *  set [Cylinders](@ref Cylinder), [Linkers](@Linker) and [MotorGhosts](@MotorGhost)
 *  with corresponding static DynamicRateChanger objects.
 */
class DRController {
    
public:
    ///Intialize all elements with corresponding DynamicRateChanger objects.
    void initialize(DynamicRateTypes& drTypes);
};


#endif
