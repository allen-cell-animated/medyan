#ifndef MEDYAN_Structure_SurfaceMesh_MembraneRegion_Hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneRegion_Hpp

#include <array>
#include <vector>
#include <memory>

#include "MathFunctions.h"
#include "Structure/Boundary.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneHierarchy.hpp"

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
    template< typename Float >
    bool contains(const mathfunc::Vec< 3, Float >& point) const {
        /**************************************************************************
        This function checks whether a certain point is in the region. If a point
        is in the region, it must satisfy all of the following.
            - within _boundary if specified
            - in any of the membranes in _hierOut
            - not in any of the membranes in _hierIn
        **************************************************************************/
        if(_boundary) {
            auto p = mathfunc::vec2Vector(point);
            if(!_boundary->within(p)) return false;
        }

        if(!_hierOut.empty()) {
            bool hierOutGood = false;
            for(auto eachHier: _hierOut) {
                if(eachHier->getMembrane()->contains(point)) {
                    hierOutGood = true;
                    break;
                }
            }
            if(!hierOutGood) return false;
        }

        for(auto eachHier: _hierIn) {
            if(eachHier->getMembrane()->contains(point)) return false;
        }

        return true;
    }

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