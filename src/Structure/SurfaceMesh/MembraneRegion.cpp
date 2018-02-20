#include "MembraneRegion.h"

#include "MembraneHierarchy.h"
#include "Membrane.h"
#include "Boundary.h"

#include "MathFunctions.h"

MembraneRegion::MembraneRegion(MembraneHierarchy* hier, bool excludeChilden):
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

bool MembraneRegion::contains(const std::array<double, 3>& point)const {
    /**************************************************************************
    This function checks whether a certain point is in the region. If a point
    is in the region, it must satisfy all of the following.
        - within _boundary if specified
        - in any of the membranes in _hierOut
        - not in any of the membranes in _hierIn
    **************************************************************************/
    if(_boundary) {
        auto p = mathfunc::array2Vector(point);
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
