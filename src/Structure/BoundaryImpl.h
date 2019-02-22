
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

#ifndef MEDYAN_BoundaryImpl_h
#define MEDYAN_BoundaryImpl_h

#include <vector>

#include "common.h"

#include "Boundary.h"

///FORWARD DECLARATIONS
class Compartment;

/// A cubic Boundary implementation.
class BoundaryCubic: public Boundary {

public:
    ///Default constructor, this will create a cube with given
    ///corners at edges of current CompartmentGrid
    BoundaryCubic(SubSystem* s, vector<BoundaryMove> move);

    virtual bool within(Compartment* C);
    virtual bool within(const vector<double>& coordinates);

    virtual double distance(const vector<double>& coordinates);

    //Qin
    virtual double lowerdistance(const vector<double>& coordinates);
    virtual double sidedistance(const vector<double>& coordinates);

    virtual double getboundaryelementcoord(int i);

    virtual void move(vector<double> dist);

    ///Returns the normal inward at this coordinate
    //rule - takes the closest wall's normal inward.
    virtual vector<double> normal(vector<double>& coordinates);

    virtual void volume();
};

/// A spherical Boundary implementation.
class BoundarySpherical: public Boundary {

public:
    ///Default constructor, will create an sphere with given diameter
    ///@param diameter - diameter of sphere
    BoundarySpherical(SubSystem* s, double diameter, vector<BoundaryMove> move);

    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<double>& coordinates);

    virtual double distance(const vector<double>& coordinates);

    //Qin
    virtual double lowerdistance(const vector<double>& coordinates);
    virtual double sidedistance(const vector<double>& coordinates);
    virtual double getboundaryelementcoord(int i){};
    ///@note - not yet implemented.
    virtual void move(vector<double> dist) {}

    ///Returns the normal inward at this coordinate
    virtual vector<double> normal(vector<double>& coordinate);

    virtual void volume();
};

/// A capsule Boundary implementation.
class BoundaryCapsule: public Boundary {

public:
    /// Default constructor, will create a capsule with given diameter, and height equal
    /// to current grid.
    /// @param diameter - diameter of capsule (will set half sphere radii as well as
    /// cylinder radius)
    BoundaryCapsule(SubSystem* s, double diameter, vector<BoundaryMove> move);

    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<double>& coordinates);

    virtual double distance(const vector<double>& coordinates);
    virtual double getboundaryelementcoord(int i){};
    //Qin
    virtual double lowerdistance(const vector<double>& coordinates);
    virtual double sidedistance(const vector<double>& coordinates);

    ///@note - Not yet implemented.
    virtual void move(vector<double> dist) {}

    ///Returns the normal inward at this coordinate
    //@note - Not yet implemented.
    virtual vector<double> normal(vector<double>& coordinate) {return vector<double>{0,0,0};}

    virtual void volume();
};

/// A cylinder Boundary implementation.
class BoundaryCylinder: public Boundary {

public:
    /// Default constructor, will create a capsule with given diameter, and height equal
    /// to current grid.
    /// @param diameter - diameter of capsule (will set half sphere radii as well as
    /// cylinder radius)
    BoundaryCylinder(SubSystem* s, double diameter, vector<BoundaryMove> move);

    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<double>& coordinates);

    virtual double distance(const vector<double>& coordinates);
    virtual double getboundaryelementcoord(int i){};
    //Qin
    virtual double lowerdistance(const vector<double>& coordinates);
    virtual double sidedistance(const vector<double>& coordinates);

    ///@note - Not yet implemented.
    virtual void move(vector<double> dist) {}

    ///Returns the normal inward at this coordinate
    //@note - Not yet implemented.
    virtual vector<double> normal(vector<double>& coordinate) {return vector<double>{0,0,0};}

    virtual void volume();
};


#endif
