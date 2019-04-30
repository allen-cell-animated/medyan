#include "Structure/SurfaceMesh/MembraneRegion.h"

#include "Structure/SurfaceMesh/MembraneHierarchy.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Boundary.h"

#include "MathFunctions.h"

MembraneRegion::MembraneRegion(MembraneHierarchy* hier, bool excludeChildren):
    _hierOut({hier})
{
    if(excludeChildren) {
        size_t n = hier->numberOfChildren();
        _hierIn.reserve(n);
        for(size_t idx = 0; idx < n; ++idx) {
            _hierIn.push_back( static_cast<MembraneHierarchy*>(hier->children(idx)) );
        }
    }
}

MembraneRegion::MembraneRegion(Boundary* b, MembraneHierarchy* parentOfExcluded):
    _boundary(b)
{
    size_t n = parentOfExcluded->numberOfChildren();
    _hierIn.reserve(n);
    for(size_t idx = 0; idx < n; ++idx) {
        _hierIn.push_back( static_cast<MembraneHierarchy*>(parentOfExcluded->children(idx)) );
    }
}

bool MembraneRegion::contains(const mathfunc::Vec3& point)const {
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

std::unique_ptr<MembraneRegion> MembraneRegion::makeByChildren(const MembraneHierarchy& hier, bool excludeChildren) {
    auto mr = unique_ptr<MembraneRegion>(new MembraneRegion());

    for(const auto& it: hier.children()) {
        auto eachHier = static_cast<MembraneHierarchy*>(it.get());
        mr->_hierOut.push_back(eachHier);
        if(excludeChildren) {
            size_t n = eachHier->numberOfChildren();
            mr->_hierIn.reserve(n);
            for(size_t idx = 0; idx < n; ++idx) {
                mr->_hierIn.push_back( static_cast<MembraneHierarchy*>(eachHier->children(idx)) );
            }
        }
    }

    return mr; // Requires copy elision
}
