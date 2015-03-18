
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

#ifndef M3SYM_BoundaryElement_h
#define M3SYM_BoundaryElement_h

#include <iostream>

#include "common.h"

#include "BoundaryElementDB.h"

#include "NeighborListDB.h"
#include "Neighbor.h"

#include "GController.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents an element of a BoundarySurface.
/*!
 * The BoundaryElement class is a representation of a BoundarySurface element, which can
 * interact with other elements in the system, including other BoundaryElements as well 
 * as [Beads] (@ref Bead) in [Filaments](@ref Filament). Together, a collection of 
 * boundary elements make up a BoundarySurface.
 */
class BoundaryElement : public Neighbor {
    
protected:
    vector<double> _coords; ///< coordinates
    
    double _kRep; ///< Repulsion constant
    double _r0; ///< Screening length
    
    
public:
    /// Default constructor
    BoundaryElement(vector<double> coords, double kRepuls, double screenLength)
        : _coords(coords), _kRep(kRepuls), _r0(screenLength) {
        
        //add to boundary element db
        BoundaryElementDB::instance()->addBoundaryElement(this);
        
        //add to neighbor list db
        NeighborListDB::instance()->addNeighbor(this);
    }
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~BoundaryElement() noexcept {
        
        //remove from boundary element db
        BoundaryElementDB::instance()->removeBoundaryElement(this);
        
        ///remove from neighbor lists
        NeighborListDB::instance()->removeNeighbor(this);
    }
    
    ///return coordinates of boundary element
    const vector<double>& getCoords() {return _coords;}
    
    /// Implement for all boundary elements
    /// Returns the distance from a given point to this boundary element
    /// @return - 1) positive number if point is within boundary element
    ///          2) Negative number if point is outside boundary element
    ///          3) Infinity if point is not in domain of this boundary element
    virtual double distance(const vector<double>& point) = 0;
    
    /// Returns stretched distance, similar to distance above
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d) = 0;
    
    /// Returns normal vector of point to plane
    virtual const vector<double> normal(const vector<double> &point) = 0;
    
    //@{
    /// Getter for parameters
    virtual double getRepulsionConst() {return _kRep;}
    virtual double getScreeningLength() {return _r0;}
    //@}
 
};

#endif
