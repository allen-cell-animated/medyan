
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Component.h"

#include "Composite.h"
#include "Visitor.h"

Composite* Component::getRoot() {
    if(isRoot())
        return static_cast<Composite*>(this);
    else
        return getParent()->getRoot();
}

bool Component::apply (Visitor &v) {
    return v.visit(this);
}