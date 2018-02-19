#include "MembraneRegion.h"

#include "MembraneHierarchy.h"
#include "Membrane.h"
#include "Boundary.h"

MembraneRegion::MembraneRegion(MembraneHierarchy* hier, bool excludeChilden):
    _hierOut({hier})
{
    if(excludeChildren) {
        // TODO
    }
}

MembraneRegion::MembraneRegion(Boundary* b, MembraneHierarchy* parentOfExcluded):
    _boundary(b)
{
    // TODO
}

bool MembraneRegion::contains(const std::array<double, 3>& point) {
    // TODO
}
