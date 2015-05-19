
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

#include "BoundaryImpl.h"

#include "BoundarySurfaceImpl.h"
#include "BoundaryElement.h"

#include "SysParams.h"

BoundaryCubic::BoundaryCubic(SubSystem* s) : Boundary(s, 3, BoundaryShape::Cube){
    
    //Get full system size (want planes to be slightly inside compartment grid)
    double zeroX = 0.1 * SysParams::Geometry().compartmentSizeX * SysParams::Geometry().NX;
    double zeroY = 0.1 * SysParams::Geometry().compartmentSizeY * SysParams::Geometry().NY;
    double zeroZ = 0.1 * SysParams::Geometry().compartmentSizeZ * SysParams::Geometry().NZ;
    
    double sysX = SysParams::Geometry().compartmentSizeX * SysParams::Geometry().NX - zeroX;
    double sysY = SysParams::Geometry().compartmentSizeY * SysParams::Geometry().NY - zeroY;
    double sysZ = SysParams::Geometry().compartmentSizeZ * SysParams::Geometry().NZ - zeroZ;
    
    //Create boundary surfaces, add to vector
    //X normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {zeroX, sysY / 2, sysZ / 2}, {1, 0, 0}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX, sysY / 2, sysZ / 2}, {-1, 0, 0}));
    
    //Y normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, zeroY, sysZ / 2}, {0, 1, 0}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY, sysZ / 2}, {0, -1, 0}));
    
    //Z normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY / 2, zeroZ}, {0, 0, 1}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY / 2, sysZ}, {0, 0, -1}));
    
}

bool BoundaryCubic::within(const vector<double>& coordinates) {
    
    // check if all planes return positive distance
    // (means in front of plane, relative to normal)
    for(auto &bs : _boundarySurfaces)
        if(bs->boundaryElements()[0]->distance(coordinates) <= 0) return false;
    return true;
}


BoundarySpherical::BoundarySpherical(SubSystem* s, double diameter)

    : Boundary(s, 3, BoundaryShape::Sphere) {
    
    double sysX = SysParams::Geometry().compartmentSizeX * SysParams::Geometry().NX;
    double sysY = SysParams::Geometry().compartmentSizeY * SysParams::Geometry().NY;
    double sysZ = SysParams::Geometry().compartmentSizeZ * SysParams::Geometry().NZ;
    
    _boundarySurfaces.emplace_back(
    new Sphere(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2));
}

bool BoundarySpherical::within(const vector<double>& coordinates) {
    
    //check if the boundary element returns a positive distance
    BoundaryElement* sphereBoundaryElement =
        _boundarySurfaces[0]->boundaryElements()[0].get();
    
    return sphereBoundaryElement->distance(coordinates) > 0;
    
}

BoundaryCapsule::BoundaryCapsule(SubSystem* s, double diameter)

    : Boundary(s, 3, BoundaryShape::Capsule) {
    
    double sysX = SysParams::Geometry().compartmentSizeX * SysParams::Geometry().NX;
    double sysY = SysParams::Geometry().compartmentSizeY * SysParams::Geometry().NY;
    double sysZ = SysParams::Geometry().compartmentSizeZ * SysParams::Geometry().NZ;

    double height = sysZ - diameter;
    
    _boundarySurfaces.emplace_back(
    new CylinderZ(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2, height));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 + height / 2}, diameter / 2, false));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 - height / 2}, diameter / 2, true));
}

bool BoundaryCapsule::within(const vector<double>& coordinates) {
    
    //check if the boundary elements return a positive distance
    for(auto &bSurface : _boundarySurfaces) {
        BoundaryElement* boundaryElement = bSurface->boundaryElements()[0].get();
        
        double dist = boundaryElement->distance(coordinates);
        if(dist <= 0) return false;
    }
    return true;
}



