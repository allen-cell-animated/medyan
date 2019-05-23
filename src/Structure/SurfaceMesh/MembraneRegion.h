#ifndef MEDYAN_MembraneRegion_h
#define MEDYAN_MembraneRegion_h

#include <array>
#include <vector>
#include <memory>

#include "MathFunctions.h"

// Forward Declarations
class MembraneHierarchy;
class Boundary;

/******************************************************************************
 * 
 * Membrane Region is a class for management of regions enclosed by membranes
 * and boundaries.
 * 
 * Note: before building the membrane region, all the related membranes must
 * have the updated geometry.
 * 
******************************************************************************/

class MembraneRegion {
private:
    Boundary* _boundary = nullptr; ///< The boundary of the playground
    std::vector<MembraneHierarchy*> _hierOut; ///< Outer limit. A point is inside the region if it is in at least 1 of the membranes.
    std::vector<MembraneHierarchy*> _hierIn;  ///< Inner limit. A point is inside the region if it is not in any of the membranes.

    /// Constructor only for internal use
    MembraneRegion() {}

public:
    /// Construct the region with one membrane
    MembraneRegion(MembraneHierarchy* hier, bool excludeChildren=false);

    /// Construct the region with the boundary
    MembraneRegion(Boundary* b): _boundary(b) {}
    MembraneRegion(Boundary* b, MembraneHierarchy* parentOfExcluded);

    /// Is point inside region
    bool contains(const mathfunc::Vec3& point)const;

    /**************************************************************************
    Getters and Setters
    **************************************************************************/
    Boundary* getBoundary()const { return _boundary; }
    const std::vector<MembraneHierarchy*>& getHierOut()const { return _hierOut; }
    const std::vector<MembraneHierarchy*>& getHierIn ()const { return _hierIn;  }

    /**************************************************************************
    Factory functions
    **************************************************************************/
    /// Create region with hier's children as outer limit
    static std::unique_ptr<MembraneRegion> makeByChildren(const MembraneHierarchy& hier, bool excludeChildren=false);
};

#endif