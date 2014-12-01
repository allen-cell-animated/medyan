
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Boundary_h
#define M3SYM_Boundary_h

#include "common.h"

#include "BoundarySurface.h"

/// BoundaryShape is a shape enumeration
enum class BoundaryShape {Cube, Capsule, Sphere};

///Boundary class is to store all boundary surfaces that are in the system
/*!
 *  The boundary class stores all [BoundarySurfaces](@ref BoundarySurface) in the given shape. Its constructors can create
 *  basic boundary shapes (for now). Eventually will be extended to more complex surfaces.
 */
class Boundary {
    
protected:
    /// Vector of boundary surfaces (could be different implementations)
    vector<unique_ptr<BoundarySurface>> _boundarySurfaces;
    BoundaryShape _shape; ///< Boundary type of this boundary
    short _nDim; ///< Dimensionality of this boundary
    
public:
    Boundary(int nDim, BoundaryShape shape) : _nDim(nDim), _shape(shape) {};
    ~Boundary() {};

    /// Get shape of boundary
    BoundaryShape getShape() {return _shape;}
    /// Get boundary surfaces
    const vector<unique_ptr<BoundarySurface>>& getBoundarySurfaces() {return _boundarySurfaces;}
    
    /// Check if coordinates are within boundary
    virtual bool within(const vector<double>& coordinates) = 0;
};


#endif
