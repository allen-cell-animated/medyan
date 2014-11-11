//
//  NeighborsManager.cpp
//  Cyto
//
//  Created by James Komianos on 11/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "NeighborsManager.h"

vector<Cylinder*> NeighborList::findNearbyCylinders(Cylinder* cylinder) {
    
    vector<Cylinder*> cylinders;
    
    ///Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
                                  SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments)
        for(auto &cyl : c->getCylinders())
            if(cylinder->getFilament() != cyl->getFilament()) {
                ///CHECK CYLINDER POSITION
                cylinders.push_back(cyl);
            }
    
    return vector<Cylinder*>(cylinders.begin(), cylinders.end());
}


void NeighborList::reset() {
    
    _list.clear(); /// clear current list
    
    //loop through all cylinders
    for(auto &cylinder : *CylinderDB::instance(CylinderDBKey())) {
        
        for(auto &nearbyCylinder : findNearbyCylinders(cylinder)) {
            
            if(cylinder->getID() < nearbyCylinder->getID()) {
                
                double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
                if(dist < _rMax && dist > _rMin) _list[cylinder].push_back(nearbyCylinder);
            }
        }
    }
}
