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

#include "CaMKIIingPoint.h"


/// A container to store a CaMKII Cylinder
/*!
 *  CaMKIICylinder class is used to manage and store the species CaMKII cylinder,
 *  off reactions, and the real position of CaMKII to be used in force field, and
 *  [NeighborLists](@ref NeighborList).
 */
class CaMKIICylinder : public Cylinder {

protected:
	CaMKIIingPoint *_camkiiPoint;
	void updateCoordinate() override;
	void updatePosition() override;


public:
	/// Constructor, initializes a cylinder
	/// composite is set to NULL and we use only one bead
	CaMKIICylinder(CaMKIIingPoint *camkiiPoint, Bead* b1, short type, int position);
	CaMKIIingPoint* getCaMKIIPointParent() {
		return _camkiiPoint;
	};
	~CaMKIICylinder();
  void addToFilamentBindingManagers();
  void removeFromFilamentBindingManagers();

};

#endif
