
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

#ifndef M3SYM_Boundary_h
#define M3SYM_Boundary_h

#include "common.h"

#include "BoundarySurface.h"

//FORWARD DECLARATIONS
class SubSystem;

/// BoundaryShape is a shape enumeration.
enum class BoundaryShape {Cube, Capsule, Sphere};

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
    
    short _nDim; ///< Dimensionality
    
public:
    Boundary(SubSystem* s, int nDim, BoundaryShape shape)
        : _subSystem(s), _shape(shape), _nDim(nDim) {};
    
    ~Boundary() {};

    /// Get shape of this boundary
    BoundaryShape getShape() {return _shape;}
    
    /// Get boundary surfaces
    const vector<unique_ptr<BoundarySurface>>& getBoundarySurfaces() {
        return _boundarySurfaces;
    }
    
    /// Check if coordinates are within boundary
    virtual bool within(const vector<double>& coordinates) = 0;
};


#endif
