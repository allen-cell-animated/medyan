
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
#include "CompartmentGrid.h"
#include "BindingManager.h"


#ifdef DYNAMICRATES
vector<Cylinder*> SubSystem::getBoundaryCylinders() {

    vector<Cylinder*> cylinders;
    
    //loop through neighbor list, construct vector
    for(auto be : BoundaryElement::getBoundaryElements()) {
        auto localCylinders = _neighborList->getNeighbors(be);
        cylinders.insert(cylinders.end(), localCylinders.begin(),
                                          localCylinders.end());
    }
    return cylinders;
}
#endif

void SubSystem::updateBindingManagers() {
    for(auto &child : _compartmentGrid->children()) {
        Compartment* c = (Compartment*)child.get();
        
        for(auto &manager : c->getFilamentBindingManagers())
            manager->updateAllPossibleBindings();
    }
}


