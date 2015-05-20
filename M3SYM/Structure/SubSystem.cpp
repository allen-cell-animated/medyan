
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

#include "SubSystem.h"

#include "BoundaryElement.h"

#ifdef DYNAMICRATES
vector<Cylinder*> SubSystem::getBoundaryCylinders() {

    vector<Cylinder*> cylinders;
    auto list = getNeighborList();
    
    //loop through neighbor list, construct vector
    for(auto be : BoundaryElement::getBoundaryElements()) {
        auto localCylinders = list->getNeighbors(be);
        cylinders.insert(cylinders.end(), localCylinders.begin(),
                                          localCylinders.end());
    }
    return cylinders;
}
#endif


