
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

protected:
	CaMKIIingPoint *_camkiiPoint;
	void updateCoordinate() override;
	void updatePosition() override;
	bool within(Cylinder* other, double dist) override;


public:
	/// Constructor, initializes a cylinder
	// composte is set to NULL and we use only one bead
	CaMKIICylinder(CaMKIIingPoint *camkiiPoint, Bead* b1, short type, int position);
	CaMKIIingPoint* getCaMKIIPointParent() {
		return _camkiiPoint;
	};
	~CaMKIICylinder();
  void addToFilamentBindingManagers();
  void removeFromFilamentBindingManagers();

};

#endif
