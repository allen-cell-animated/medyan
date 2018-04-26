
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

#include "BoundaryImpl.h"

#include "BoundarySurfaceImpl.h"
#include "BoundaryElement.h"

#include "Compartment.h"

#include "core/controller/GController.h"
#include "SysParams.h"

BoundaryCubic::BoundaryCubic(SubSystem* s, BoundaryMove move)

    : Boundary(s, 3, BoundaryShape::Cube, move){
    
    //Get full system size (want planes to be slightly inside compartment grid)
    double zeroX = 25;
    double zeroY = 25;
    double zeroZ = 25;
    
    double sysX = GController::getSize()[0] - zeroX;
    double sysY = GController::getSize()[1] - zeroY;
    double sysZ = GController::getSize()[2] - zeroZ;
    
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

bool BoundaryCubic::within(Compartment* C) {
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        //go half a compartment dist in the direction of normal
        auto coordinate = C->coordinates();
        
        //initial check of coord
        if(be->distance(coordinate) > 0) continue;
        
        //if not, see if any part is in bounds
        auto normal = be->normal(coordinate);
        
        auto effCoordinate =
        {coordinate[0] + SysParams::Geometry().compartmentSizeX * normal[0] / 2,
         coordinate[1] + SysParams::Geometry().compartmentSizeY * normal[1] / 2,
         coordinate[2] + SysParams::Geometry().compartmentSizeZ * normal[2] / 2};
        
        //check if this is within boundary
        if(be->distance(effCoordinate) <= 0)
            return false;
    
    }
    return true;
}


bool BoundaryCubic::within(const vector<double>& coordinates) {
    
    // check if all planes return positive distance
    // (means in front of plane, relative to normal)
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        if(be->distance(coordinates) <= 0)
            return false;
    }
    return true;
}

double BoundaryCubic::distance(const vector<double>& coordinates) {
    
    // loop through, get smallest non-negative distance
    double smallestDist = numeric_limits<double>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->distance(coordinates);
        
        if(dist < 0) continue;
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

vector<double> BoundaryCubic::normal(vector<double>& coordinates) {
    
    // loop through, get smallest non-negative distance
    BoundaryElement* closestPlane = nullptr;
    double smallestDist = numeric_limits<double>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->distance(coordinates);
        
        if(dist < 0) continue;
        if(dist < smallestDist) {
            smallestDist = dist;
            closestPlane = be;
        }
        
    }
    //return normal of plane
    return closestPlane->normal(coordinates);
}

void BoundaryCubic::move(double dist) {
    
    //do nothing
    if(_move == BoundaryMove::None) return;
    
    //move the top plane
    else if(_move == BoundaryMove::Top) {
        
        for(auto &bs : _boundarySurfaces) {
            
            auto be = bs->boundaryElements()[0].get();
            
            if(be->normal({0,0,0})[2] < 0) {
                
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] + dist});
                return;
            }
        }
    }
    
    else if(_move == BoundaryMove::All) {
        
        for(auto &bs : _boundarySurfaces) {
            
            auto be = bs->boundaryElements()[0].get();
            
            
            if(be->normal({0,0,0})[0] > 0)
                be->updateCoords({be->_coords[0] - dist, be->_coords[1], be->_coords[2]});
                
            else if (be->normal({0,0,0})[0] < 0)
                be->updateCoords({be->_coords[0] + dist, be->_coords[1], be->_coords[2]});
                
            else if (be->normal({0,0,0})[1] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1] - dist, be->_coords[2]});
                
            else if (be->normal({0,0,0})[1] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1] + dist, be->_coords[2]});
                
            else if (be->normal({0,0,0})[2] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] - dist});
                
            else if (be->normal({0,0,0})[2] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] + dist});

        }
    }
}


BoundarySpherical::BoundarySpherical(SubSystem* s, double diameter, BoundaryMove move)

    : Boundary(s, 3, BoundaryShape::Sphere, move) {
    
    double sysX = GController::getSize()[0];
    double sysY = GController::getSize()[1];
    double sysZ = GController::getSize()[2];
        
    _boundarySurfaces.emplace_back(
    new Sphere(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2));
}

bool BoundarySpherical::within(Compartment* C) {
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        //go half a compartment dist in the direction of normal
        auto coordinate = C->coordinates();
        
        //initial check of coord
        if(be->distance(coordinate) > 0) continue;
        
        //if not, see if any part is in bounds
        auto normal = be->normal(coordinate);
        
        auto effCoordinate =
        {coordinate[0] + SysParams::Geometry().compartmentSizeX * normal[0] / 2,
         coordinate[1] + SysParams::Geometry().compartmentSizeY * normal[1] / 2,
         coordinate[2] + SysParams::Geometry().compartmentSizeZ * normal[2] / 2};
        
        //check if this is within boundary
        if(be->distance(effCoordinate) <= 0)
            return false;
        
    }
    return true;
}

bool BoundarySpherical::within(const vector<double>& coordinates) {
    
    //check if the boundary element returns a positive distance
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    return be->distance(coordinates) > 0;
    
}

double BoundarySpherical::distance(const vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    double dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<double>::infinity();
}


vector<double> BoundarySpherical::normal(vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    return be->normal(coordinates);
}


BoundaryCapsule::BoundaryCapsule(SubSystem* s, double diameter, BoundaryMove move)

    : Boundary(s, 3, BoundaryShape::Capsule, move) {
    
    double sysX = GController::getSize()[0];
    double sysY = GController::getSize()[1];
    double sysZ = GController::getSize()[2];

    double height = sysZ - diameter;
    
    _boundarySurfaces.emplace_back(
    new CylinderZ(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2, height));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 + height / 2}, diameter / 2, false));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 - height / 2}, diameter / 2, true));
}

bool BoundaryCapsule::within(Compartment* C) {
    
    //just calls regular within for now
    return within(C->coordinates());
}


bool BoundaryCapsule::within(const vector<double>& coordinates) {
    
    //check if the boundary elements return a positive distance
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->distance(coordinates);
        if(dist <= 0) return false;
    }
    return true;
}

double BoundaryCapsule::distance(const vector<double>& coordinates) {
    
    // loop through, get smallest non-negative distance
    double smallestDist = numeric_limits<double>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->distance(coordinates);
        
        if(dist < 0) continue;
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

