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
    
    ///Create points of cube
    std::vector<double> p0 = {zeroX, zeroY, zeroZ};
    std::vector<double> p1 = {sysX, zeroY, zeroZ};
    std::vector<double> p2 = {sysX, sysY, zeroZ};
    std::vector<double> p3 = {zeroX, sysY, zeroZ};
    
    std::vector<double> p4 = {zeroX, zeroY, sysZ};
    std::vector<double> p5 = {sysX, zeroY, sysZ};
    std::vector<double> p6 = {sysX, sysY, sysZ};
    std::vector<double> p7 = {zeroX, sysY, sysZ};
    
    std::vector<int> numDivisions = {SystemParameters::Geometry().compartmentSizeX * 10,
                                     SystemParameters::Geometry().compartmentSizeY * 10}
    
    ///Create boundary surfaces, add to vector
    ///Y normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p0, p1, p5, p4}, numDivisions, 1)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p3, p2, p6, p7}, numDivisions, 1)));
    
    ///X normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p0, p3, p7, p4}, numDivisions, 0)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p1, p2, p6, p5}, numDivisions, 0)));
    
    ///Z normal planes
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p0, p1, p2, p3}, numDivisions, 2)));
    _boundarySurfaces.push_back(make_unique<BoundarySurface>
        (new BasicPlane({p4, p5, p6, p7}, numDivisions, 2)));
    
}