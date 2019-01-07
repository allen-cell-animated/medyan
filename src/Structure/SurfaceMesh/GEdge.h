#ifndef MEDYAN_GEdge_h
#define MEDYAN_GEdge_h

#include "MathFunctions.h"

struct GEdge {

    double length; // Current length
    double sLength; // Temporarily store the stretched length

    mathfunc::Vec3 pseudoUnitNormal; // The pseudo unit normal vector at the edge pointing outward.
    mathfunc::Vec3 sPseudoUnitNormal; // The pseudo normal under stretch

    // Auxilliary getters
    template< bool stretched > auto& getLength() { return stretched ? sLength : length; }
    template< bool stretched > auto& getPseudoUnitNormal() { return stretched ? sPseudoUnitNormal : pseudoUnitNormal; }
    
};


#endif
