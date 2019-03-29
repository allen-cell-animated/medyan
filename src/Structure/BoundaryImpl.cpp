
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

BoundaryCubic::BoundaryCubic(SubSystem* s, vector<BoundaryMove> move)

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

    volume();
    
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

//Qin
double BoundaryCubic::lowerdistance(const vector<double>& coordinates) {
    
    // loop through, get smallest non-negative distance
    double smallestDist = numeric_limits<double>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->lowerdistance(coordinates);
        
        if(dist < 0) continue;
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

//Qin, the same as lowerdistance for now
double BoundaryCubic::sidedistance(const vector<double>& coordinates) {
    
    // loop through, get smallest non-negative distance
    double smallestDist = numeric_limits<double>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        double dist = be->lowerdistance(coordinates);
        
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

double BoundaryCubic::getboundaryelementcoord(int bidx) {

    for(auto &bs : _boundarySurfaces) {

        auto be = bs->boundaryElements()[0].get();

        if(be->normal({0,0,0})[0] > 0 && bidx == 0)
            return be->_coords[0];

        else if (be->normal({0,0,0})[0] < 0 && bidx == 1)
            return be->_coords[0];

        else if (be->normal({0,0,0})[1] > 0 && bidx == 2)
            return be->_coords[1];

        else if (be->normal({0,0,0})[1] < 0 && bidx == 3)
            return be->_coords[1];

        else if (be->normal({0,0,0})[2] > 0&& bidx == 4)
            return be->_coords[2];

        else if (be->normal({0,0,0})[2] < 0&&  bidx == 5 )
            return be->_coords[2];

    }

};

void BoundaryCubic::move(vector<double> dist) {

    //do nothing
    if(_move.size() ==0 ) return;
    else if(_move.size() == 1 && _move[0] == BoundaryMove::All) {

        for(auto &bs : _boundarySurfaces) {

            auto be = bs->boundaryElements()[0].get();

            //right
            if(be->normal({0,0,0})[0] < 0)
                be->updateCoords({be->_coords[0] + dist[1], be->_coords[1],
                                  be->_coords[2]});
                //left
            else if (be->normal({0,0,0})[0] > 0)
                be->updateCoords({be->_coords[0] + dist[0], be->_coords[1],
                                  be->_coords[2]});
                //top
            else if (be->normal({0,0,0})[1] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1] + dist[3],
                                  be->_coords[2]});
                //bottom
            else if (be->normal({0,0,0})[1] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1] + dist[2],
                                  be->_coords[2]});
                //front
            else if (be->normal({0,0,0})[2] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] +
                                                                  dist[5]});
                //back
            else if (be->normal({0,0,0})[2] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] +
                                                                  dist[4]});

        }
    }
    else {
        for (auto bm:_move) {
            if (bm == BoundaryMove::None) return;
                //move the left plane
            else if (bm == BoundaryMove::Left) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[0] > 0) {

                        be->updateCoords(
                                {be->_coords[0] + dist[0], be->_coords[1], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Right) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[0] < 0) {

                        be->updateCoords(
                                {be->_coords[0] + dist[1], be->_coords[1], be->_coords[2] });
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Front) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[1] > 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1] + dist[2], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Back) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[1] < 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1] + dist[3], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Bottom) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[2] > 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1], be->_coords[2] + dist[4]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Top) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[2] < 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1], be->_coords[2] + dist[5]});
                        return;
                    }
                }
            }
        }
    }
}

void BoundaryCubic::volume(){
    // array
    double spanarray[3]={0.0,0.0,0.0}; //stores span along x y and z axes.
    for(auto &bs : _boundarySurfaces) {

        auto be = bs->boundaryElements()[0].get();

        auto coord = be->_coords;

        //right
        if(be->normal({0,0,0})[0] > 0)
            spanarray[0] -= coord[0];
        //left
        else if (be->normal({0,0,0})[0] < 0)
            spanarray[0] += coord[0];
        //top
        else if (be->normal({0,0,0})[1] > 0)
            spanarray[1] -= coord[1];
        //bottom
        else if (be->normal({0,0,0})[1] < 0)
            spanarray[1] += coord[1];
        //front
        else if (be->normal({0,0,0})[2] > 0)
            spanarray[2] -= coord[2];
        //back
        else if (be->normal({0,0,0})[2] < 0)
            spanarray[2] += coord[2];

    }

    Boundary::systemvolume = (50 + spanarray[0]) * (50 + spanarray[1]) * (50 +
            spanarray[2]);
}


BoundarySpherical::BoundarySpherical(SubSystem* s, double diameter, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Sphere, move) {
    
    double sysX = GController::getSize()[0];
    double sysY = GController::getSize()[1];
    double sysZ = GController::getSize()[2];
        
    _boundarySurfaces.emplace_back(
    new Sphere(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2));
    volume();
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

//Qin, the same as distance
double BoundarySpherical::lowerdistance(const vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    double dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<double>::infinity();
}

double BoundarySpherical::sidedistance(const vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    double dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<double>::infinity();
}


vector<double> BoundarySpherical::normal(vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    return be->normal(coordinates);
}

void BoundarySpherical::volume(){
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    double r[1];
    be->elementeqn(r);
    Boundary::systemvolume = 4/3 * 22/7 * r[0] * r[0] * r[0];
}


BoundaryCapsule::BoundaryCapsule(SubSystem* s, double diameter, vector<BoundaryMove> move)

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
    volume();
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

//Qin, the same as distance
double BoundaryCapsule::lowerdistance(const vector<double>& coordinates) {
    
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

double BoundaryCapsule::sidedistance(const vector<double>& coordinates) {
    
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

void BoundaryCapsule::volume(){
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    double dim[2];
    be->elementeqn(dim);
    double r = dim[0];
    double h = dim[1];
    double pi = 22/7;
    Boundary::systemvolume = (h + 4/3 * r) * pi * r * r;
}

///Qin------------------


BoundaryCylinder::BoundaryCylinder(SubSystem* s, double diameter, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Cylinder, move) {
    
    double sysX = GController::getSize()[0];
    double sysY = GController::getSize()[1];
    double sysZ = GController::getSize()[2];
    
    double height = sysZ;
    
    _boundarySurfaces.emplace_back(
                                   new CylinderXYZ(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2, height));
    volume();

}

bool BoundaryCylinder::within(Compartment* C) {

    // project compartment to a 2D cylinderical coordinate
    double comX = GController::getCompartmentSize()[0];
    double comY = GController::getCompartmentSize()[1];
    auto r = SysParams::Boundaries().diameter / 2;
    auto x = C->coordinates()[0] - r;
    auto y = C->coordinates()[1] - r;

    auto r1 = sqrt((x - comX / 2) * (x - comX / 2) + (y - comY / 2) * (y - comY / 2));
    auto r2 = sqrt((x + comX / 2) * (x + comX / 2) + (y - comY / 2) * (y - comY / 2));
    auto r3 = sqrt((x - comX / 2) * (x - comX / 2) + (y + comY / 2) * (y + comY / 2));
    auto r4 = sqrt((x + comX / 2) * (x + comX / 2) + (y + comY / 2) * (y + comY / 2));
    
    //cout << "x= " << C->coordinates()[0] << " y = " << C->coordinates()[1] << endl;
    
    if (r1 < r || r2 < r || r3 < r || r4 < r) return true;
    else return false;
    
    //return within(C->coordinates());
}


bool BoundaryCylinder::within(const vector<double>& coordinates) {
    
    //check if the boundary elements return a positive distance
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();

        double dist = be->distance(coordinates);
        if(dist <= 0) return false;
    }
    return true;
}

double BoundaryCylinder::distance(const vector<double>& coordinates) {
    
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

//Qin
//lower distance should only return the distance between beads and the lower boundary
double BoundaryCylinder::lowerdistance(const vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    double dist = be->lowerdistance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<double>::infinity();
}

double BoundaryCylinder::sidedistance(const vector<double>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    double dist = be->sidedistance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<double>::infinity();
}

void BoundaryCylinder::volume(){
    double radius = GController::getSize()[0];
    double sysZ = GController::getSize()[2];

    Boundary::systemvolume = radius * radius * 3.14159 * sysZ;
}

