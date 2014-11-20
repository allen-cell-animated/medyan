//
//  BoundaryElement.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElement__
#define __Cyto__BoundaryElement__

#include <vector>
#include <iostream>

#include "common.h"
#include "Gcontroller.h"

///FORWARD DECLARATIONS
class Bead;

///BoundaryElement class represents an element of a boundary surface
/*!
 * The BoundaryElement class is a representation of a boundary surface element, which can interact
 * with other elements in the system, including other BoundaryElements as well as Beads in filaments.
 * Together, a collection of boundary elements make up a BoundarySurface.
 */
class BoundaryElement {
    
protected:
    vector<Bead*> _beads; ///< Beads that this boundary element could interact with
    vector<double> _coords; ///< coordinates of this boundary element
    
public:
    ///Default constructor
    BoundaryElement(vector<double> coords) : _coords(coords) {}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~BoundaryElement() noexcept {}
    
    ///add a bead to list of interacting beads
    void addBead(Bead* b) {_beads.push_back(b);}
    ///Remove a bead from list of interacting beads
    ///@note does nothing if bead is not in interacting list already
    void removeBead(Bead* b) {
        auto it = find(_beads.begin(), _beads.end(), b);
        if(it != _beads.end()) _beads.erase(it);
    }
    
    ///return coordinates of boundary element
    const vector<double>& getCoords() {return _coords;}
    ///return set of beads
    const vector<Bead*>& getBeads() {return _beads;}
    
    ///Implement for all boundary elements
    ///Returns the distance from a given point to this boundary element
    ///@return - 1) positive number if point is within boundary element
    ///          2) Negative number if point is outside boundary element
    ///          3) Infinity if point is not in domain of this boundary element
    virtual double distance(const vector<double>& point) = 0;
    
    ///Returns stretched distance, similar to distance above
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d) = 0;
    
    ///Returns normal vector of point to plane
    virtual const vector<double> normal(const vector<double> &point) = 0;
    
    ///Getters for parameters
    virtual double getRepulsionConst() = 0;
    virtual double getScreeningLength() = 0;
 
};



#endif /* defined(__Cyto__BoundaryElement__) */
