
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Boundary_h
#define MEDYAN_Boundary_h

#include "common.h"

#include "BoundarySurface.h"

//FORWARD DECLARATIONS
class SubSystem;

/// BoundaryShape is a shape enumeration.
enum class BoundaryShape {Cube, Capsule, Sphere};

/// BoundaryMove is a enum describing the movement of a boundary.
enum class BoundaryMove {None, Top, All};

/// To store all [BoundarySurfaces](@ref BoundarySurface) that are in the SubSystem.
/*!
 *  The boundary class stores all [BoundarySurfaces](@ref BoundarySurface) in the given 
 *  shape. Its constructors can create basic boundary shapes (for now). Eventually will 
 *  be extended to more complex surfaces.
 */
class Boundary {
    
protected:
    SubSystem* _subSystem; ///< SubSystem ptr
    
    /// Vector of boundarysurfaces (could be different implementations)
    vector<unique_ptr<BoundarySurface>> _boundarySurfaces;
    
    BoundaryShape _shape; ///< Shape of boundary
    BoundaryMove _move;   ///< Movement of boundary
    
    short _nDim; ///< Dimensionality
    
public:
    Boundary(SubSystem* s, int nDim, BoundaryShape shape, BoundaryMove move)
        : _subSystem(s), _shape(shape), _move(move), _nDim(nDim) {};
    
    ~Boundary() {};

    /// Get shape of this boundary
    BoundaryShape getShape() {return _shape;}
    
    /// Get boundary surfaces
    const vector<unique_ptr<BoundarySurface>>& getBoundarySurfaces() {
        return _boundarySurfaces;
    }
    
    /// Check if coordinates are within boundary
    virtual bool within(const vector<double>& coordinates) = 0;
    
    /// Check if a compartment is within boundary
    /// @note - this checks if ANY part of the compartment volume
    ///         is within the boundary.
    virtual bool within(Compartment* C) = 0;
    
    /// Get the distance from the boundary. Returns the distance from
    /// closest boundary element in the boundary.
    /// Will return infinity if outside of the boundary.
    virtual double distance(const vector<double>& coordinates) = 0;
    
    ///Move a given part of a boundary a given distance
    ///@note a negative distance denotes movement towards the center of the grid.
    virtual void move(double dist) = 0;
    
    //Give a normal to the plane (pointing inward) at a given point
    virtual vector<double> normal(vector<double>& coordinates) = 0;
    
};


#endif
