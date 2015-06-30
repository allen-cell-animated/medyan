
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

#ifndef M3SYM_BoundaryImpl_h
#define M3SYM_BoundaryImpl_h

#include <vector>

#include "common.h"

#include "Boundary.h"

/// A cubic Boundary implementation.
class BoundaryCubic: public Boundary {
    
public:
    ///Default constructor, this will create a cube with given corners at edges of
    ///current CompartmentGrid
    BoundaryCubic(SubSystem* s);
    
    virtual bool within(const vector<double>& coordinates);
    virtual double distance(const vector<double>& coordinates);
};

/// A spherical Boundary implementation.
class BoundarySpherical: public Boundary {
    
public:
    ///Default constructor, will create an sphere with given diameter
    ///@param diameter - diameter of sphere
    BoundarySpherical(SubSystem* s, double diameter);
    
    virtual bool within(const vector<double>& coordinates);
    virtual double distance(const vector<double>& coordinates);
};

/// A capsule Boundary implementation.
class BoundaryCapsule: public Boundary {
    
public:
    /// Default constructor, will create a capsule with given diameter, and height equal
    /// to current grid.
    /// @param diameter - diameter of capsule (will set half sphere radii as well as
    /// cylinder radius)
    BoundaryCapsule(SubSystem* s, double diameter);
    
    virtual bool within(const vector<double>& coordinates);
    virtual double distance(const vector<double>& coordinates);
};


#endif
