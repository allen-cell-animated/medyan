//
//  Boundary.h
//  Cyto
//
//  Created by James Komianos on 8/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Boundary__
#define __Cyto__Boundary__

#include <iostream>

#include "common.h"
#include "BoundarySurface.h"

///Boundary type enumeration
enum class BoundaryShape {Cube, Capsule, Sphere};

///Boundary class is to store all boundary surfaces that are in the system
/*!
 *  The boundary class stores all boundary surfaces in the given shape. Its constructors can create
 *  basic boundary shapes (for now). Eventually will be extended to more complex surfaces.
 */
class Boundary {
    
protected:
    vector<unique_ptr<BoundarySurface>> _boundarySurfaces;
    ///<vector of boundary surfaces (could be different implementations)
    BoundaryShape _shape; ///<Boundary type of this boundary
    short _nDim; ///<Dimensionality of this boundary
    
public:
    ///not much for now
    
    ///Constructor and destructor
    Boundary(int nDim, BoundaryShape shape) : _nDim(nDim), _shape(shape) {};
    ~Boundary() {};

    //get shape
    BoundaryShape getShape() {return _shape;}
    ///get boundary surfaces
    const vector<unique_ptr<BoundarySurface>>& getBoundarySurfaces() {return _boundarySurfaces;}
    
    ///Check if coordinates are within boundary 
    virtual bool within(const vector<double> coordinates) = 0;
};


#endif /* defined(__Cyto__Boundary__) */
