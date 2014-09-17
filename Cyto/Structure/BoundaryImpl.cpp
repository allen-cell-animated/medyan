//
//  BoundaryImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryImpl.h"

BoundaryCubic::BoundaryCubic() : Boundary(3, BoundaryShape::Cube){
    
    ///Get full system size (want planes to be slightly inside compartment grid)
    
    double zeroX = 0.01 * SystemParameters::Geometry().compartmentSizeX;
    double zeroY = 0.01 * SystemParameters::Geometry().compartmentSizeY;
    double zeroZ = 0.01 * SystemParameters::Geometry().compartmentSizeZ;
    
    double sysX = SystemParameters::Geometry().compartmentSizeX * SystemParameters::Geometry().NX - zeroX;
    double sysY = SystemParameters::Geometry().compartmentSizeY * SystemParameters::Geometry().NY - zeroY;
    double sysZ = SystemParameters::Geometry().compartmentSizeZ * SystemParameters::Geometry().NZ - zeroZ;
    
    ///Create boundary surfaces, add to vector
    ///X normal planes
    _boundarySurfaces.emplace_back(new Plane({0, sysY / 2, sysZ / 2}, {1, 0, 0}));
    _boundarySurfaces.emplace_back(new Plane({sysX, sysY / 2, sysZ / 2}, {-1, 0, 0}));
    
    ///Y normal planes
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, 0, sysZ / 2}, {0, 1, 0}));
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY, sysZ / 2}, {0, -1, 0}));
    
    ///Z normal planes
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY / 2, 0}, {0, 0, 1}));
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY / 2, sysZ}, {0, 0, -1}));
    
}