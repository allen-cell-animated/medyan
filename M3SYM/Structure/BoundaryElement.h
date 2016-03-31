
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryElement_h
#define M3SYM_BoundaryElement_h

#include <iostream>

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Neighbor.h"
#include "Component.h"

#include "GController.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents an element of a BoundarySurface.
/*!
 * The BoundaryElement class is a representation of a BoundarySurface element, which can
 * interact with other elements in the system, including other BoundaryElements as well 
 * as [Beads] (@ref Bead) in [Filaments](@ref Filament) and [Bubbles](@ref Bubble). 
 * Together, a collection of boundary elements make up a BoundarySurface. 
 *
 * Extending the Neighbor class, all instances can be kept in 
 * [NeighborLists](@ref NeighborList).
 */
class BoundaryElement : public Component, public Trackable, public Neighbor {

friend class BoundaryCubic;
friend class BoundarySpherical;
friend class BoundaryCapsule;
    
private:
    static Database<BoundaryElement*> _boundaryElements;
    ///< Collection of boundary elements in SubSystem
    
protected:
    
    vector<double> _coords; ///< coordinates
    
    double _kRep; ///< Repulsion constant
    double _r0; ///< Screening length
    
public:
    /// Default constructor
    BoundaryElement(vector<double> coords, double kRepuls, double screenLength)
    
        : Trackable(false, false, false, true),
          _coords(coords), _kRep(kRepuls), _r0(screenLength) {}
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~BoundaryElement() noexcept {}
    
    ///return coordinates of boundary element
    const vector<double>& getCoords() {return _coords;}
    
    ///update the coordinates of the boundary element
    virtual void updateCoords(const vector<double> newCoords) = 0;
    
    /// Implement for all boundary elements
    /// Returns the distance from a given point to this boundary element
    /// @return - 1) positive number if point is within boundary element
    ///           2) Negative number if point is outside boundary element
    ///           3) Infinity if point is not in domain of this boundary element
    virtual double distance(const vector<double>& point) = 0;
    
    /// Returns stretched distance, similar to distance above
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d) = 0;
    
    /// Returns normal vector of point to plane
    virtual const vector<double> normal(const vector<double> &point) = 0;
    
    //@{
    /// Getter for mechanical parameters
    virtual double getRepulsionConst() {return _kRep;}
    virtual double getScreeningLength() {return _r0;}
    //@}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _boundaryElements.addElement(this);}
    virtual void removeFromSubSystem() {_boundaryElements.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<BoundaryElement*>& getBoundaryElements() {
        return _boundaryElements.getElements();
    }
    /// Get the number of boundary elements in this system
    static int numBoundaryElements() {
        return _boundaryElements.countElements();
    }
    
    virtual void printSelf();
};

#endif
