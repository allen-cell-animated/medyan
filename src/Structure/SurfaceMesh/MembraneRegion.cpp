#include "Structure/SurfaceMesh/MembraneRegion.h"

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

std::unique_ptr<MembraneRegion> MembraneRegion::makeByChildren(const MembraneHierarchy& hier, bool excludeChildren) {
    auto mr = unique_ptr<MembraneRegion>(new MembraneRegion());

    for(const auto& it: hier.children()) {
        auto eachHier = static_cast<MembraneHierarchy*>(it.get());
        if(eachHier->getMembrane()->isClosed()) {
            mr->_hierOut.push_back(eachHier);

            if(excludeChildren) {
                size_t n = eachHier->numberOfChildren();
                for(size_t idx = 0; idx < n; ++idx) {
                    auto childHier = static_cast< MembraneHierarchy* >(eachHier->children(idx));
                    if(childHier->getMembrane()->isClosed())
                        mr->_hierIn.push_back( childHier );
                }
            }
        }
    }

    return mr; // Requires copy elision
}
