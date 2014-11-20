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

#include "NeighborListDB.h"
#include "Neighbor.h"
#include "Gcontroller.h"

///FORWARD DECLARATIONS
class Bead;

///BoundaryElement class represents an element of a boundary surface
/*!
 * The BoundaryElement class is a representation of a boundary surface element, which can interact
 * with other elements in the system, including other BoundaryElements as well as Beads in filaments.
 * Together, a collection of boundary elements make up a BoundarySurface.
 */
class BoundaryElement : public Neighbor {
    
protected:
    vector<double> _coords; ///< coordinates of this boundary element
    
public:
    ///Default constructor
    BoundaryElement(vector<double> coords) : _coords(coords) {
        ///add to neighbor list db
        NeighborListDB::instance()->addNeighbor(this);
    }
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~BoundaryElement() noexcept {}
    
    ///return coordinates of boundary element
    const vector<double>& getCoords() {return _coords;}
    
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
