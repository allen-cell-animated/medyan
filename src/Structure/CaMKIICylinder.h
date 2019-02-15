
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

#ifndef MEDYAN_CaMKII_Cylinder_h
#define MEDYAN_CaMKII_Cylinder_h

#include <iostream>
#include "Cylinder.h"

//FORWARD DECLARATIONS
class CaMKIIingPoint;

/// A container to store a CaMKII Cylinder
/*!
 *  TODO CaMKII update doc:
 *  Cylinder class is used to manage and store a MCylinder and CCylinder.
 *  Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with 
 *  all instances can be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be 
 *  kept in [NeighborLists](@ref NeighborList).
 */
class CaMKIICylinder : public Cylinder {

private:
	CaMKIIingPoint *camkiiPoint;

public:
    /// Constructor, initializes a cylinder
	CaMKIICylinder(CaMKIIingPoint *camkiiPoint, Composite* parent, Bead* b1, Bead* b2, short type, int position,
             bool extensionFront = false, bool extensionBack  = false,
			 bool initialization = false):Cylinder(parent, b1, b2, type, position, extensionFront, extensionBack, initialization), camkiiPoint(camkiiPoint){
	}

};

#endif
