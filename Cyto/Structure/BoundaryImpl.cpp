//
//  BoundaryImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryImpl.h"

#include "BoundarySurfaceImpl.h"
#include "SystemParameters.h"

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
    _boundarySurfaces.emplace_back(new Plane({zeroX, sysY / 2, sysZ / 2}, {1, 0, 0}));
    _boundarySurfaces.emplace_back(new Plane({sysX, sysY / 2, sysZ / 2}, {-1, 0, 0}));
    
    ///Y normal planes
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, zeroY, sysZ / 2}, {0, 1, 0}));
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY, sysZ / 2}, {0, -1, 0}));
    
    ///Z normal planes
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY / 2, zeroZ}, {0, 0, 1}));
    _boundarySurfaces.emplace_back(new Plane({sysX / 2, sysY / 2, sysZ}, {0, 0, -1}));
    
}

bool BoundaryCubic::within(const vector<double> coordinates) {
    
    ///check if all planes return positive distance (means in front of plane, relative to normal)
    for(auto &bs : _boundarySurfaces)
        if(bs->boundaryElements()[0]->distance(coordinates) <= 0) return false;
    return true;
}


BoundarySpherical::BoundarySpherical() : Boundary(3, BoundaryShape::Sphere) {
    
    
    double sysX = SystemParameters::Geometry().compartmentSizeX * SystemParameters::Geometry().NX;
    double sysY = SystemParameters::Geometry().compartmentSizeY * SystemParameters::Geometry().NY;
    double sysZ = SystemParameters::Geometry().compartmentSizeZ * SystemParameters::Geometry().NZ;
    
    _boundarySurfaces.emplace_back(new Sphere({sysX / 2, sysY / 2, sysZ / 2}, sysX / 2));
}

bool BoundarySpherical::within(const vector<double> coordinates) {
    
    BoundaryElement* sphereBoundaryElement = _boundarySurfaces[0]->boundaryElements()[0].get();
    return sphereBoundaryElement->distance(coordinates) > 0;
    
}

