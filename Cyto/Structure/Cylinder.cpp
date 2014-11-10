//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
#include <exception>
#include "Cylinder.h"

#include "ChemManager.h"
#include "Composite.h"
#include "GController.h"
#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

Cylinder::Cylinder(Filament* f, Bead* b1, Bead* b2, bool extensionFront, bool extensionBack, bool creation) {
    
    ///Set beads
    _b1 = b1;
    _b2 = b2;
    setFilament(f);
    
    ///check if were still in same compartment
    coordinate = MidPointCoordinate(b1->coordinate, b2->coordinate, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    vector<Cylinder*> neighbors;
    ///If not initaliziation, get a neighbors list
    if(creation || extensionFront || extensionBack)  neighbors = findNearbyCylinders();
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(
        ChemManager::createCCylinder(ChemManagerCylinderKey(), f, _compartment, extensionFront, extensionBack, creation));
    _cCylinder->setCylinder(this);
    
    
    if(creation || extensionFront || extensionBack) {
        vector<CCylinder*> cNeighbors(neighbors.size());
        transform(neighbors.begin(),neighbors.end(),cNeighbors.begin(), [](Cylinder *c){return c->getCCylinder();});
        
        ///Update filament reactions, only if not initialization
        ChemManager::updateCCylinder(ChemManagerCylinderKey(), _cCylinder.get(), cNeighbors);
    }
    
#endif

#ifdef MECHANICS
    double eqLength;
    
    ///set eqLength according to cylinder size
    if(extensionFront || extensionBack) eqLength = SystemParameters::Geometry().monomerSize;
    else if(creation) eqLength = SystemParameters::Geometry().monomerSize * 2;
    else eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = unique_ptr<MCylinder>(new MCylinder(eqLength));
    
    _mCylinder->setCylinder(this);
    _mCylinder->setCoordinate(coordinate);
    
    ///Update neighbors list
    if(creation || extensionFront || extensionBack) {
        vector<MCylinder*> mNeighbors(neighbors.size());
        transform(neighbors.begin(),neighbors.end(),mNeighbors.begin(), [](Cylinder *c){return c->getMCylinder();});
        
        _mCylinder->updateExVolNeighborsList(mNeighbors);
    }
#endif
    
    ///add to compartment
    _compartment->addCylinder(this);
}

Cylinder::~Cylinder() {

    ///remove from compartment
    _compartment->removeCylinder(this);
}

vector<Cylinder*> Cylinder::findNearbyCylinders() {
    
    vector<Cylinder*> cylinders;
    
    ///Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    GController::findCompartments(coordinate, _compartment,
        SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments)
        for(auto &cyl : c->getCylinders())
            if(_pFilament != cyl->getFilament()) cylinders.push_back(cyl);

    return vector<Cylinder*>(cylinders.begin(), cylinders.end());
}


void Cylinder::updatePosition() {

    ///check if were still in same compartment, set new position
    coordinate = MidPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        ///remove from old compartment, add to new
        _compartment->removeCylinder(this);
        _compartment = c;
        _compartment->addCylinder(this);
        
#ifdef CHEMISTRY
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
#endif
    }
    ///Get a neighbors list
    auto neighbors = findNearbyCylinders();
    
    ///Update filament reactions
#ifdef CHEMISTRY
    vector<CCylinder*> cNeighbors(neighbors.size());
    transform(neighbors.begin(),neighbors.end(),cNeighbors.begin(), [](Cylinder *c){return c->getCCylinder();});
    
    ChemManager::updateCCylinder(ChemManagerCylinderKey(), _cCylinder.get(), cNeighbors);
#endif
    
    ///update exvol neighbors list
#ifdef MECHANICS
    _mCylinder->setCoordinate(coordinate);
    
    vector<MCylinder*> mNeighbors(neighbors.size());
    transform(neighbors.begin(),neighbors.end(),mNeighbors.begin(), [](Cylinder *c){return c->getMCylinder();});
    
    _mCylinder->updateExVolNeighborsList(mNeighbors);
#endif
    
}


