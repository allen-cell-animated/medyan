//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
#include <exception>
#include "Cylinder.h"

#include "ChemInitializer.h"
#include "Composite.h"
#include "GController.h"
#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

Cylinder::Cylinder(Filament* pf, Bead* firstBead, Bead* secondBead, bool extensionFront, bool extensionBack, bool creation) {
    
    ///Set beads
    _pFirst = firstBead;
    _pSecond = secondBead;
    setFilament(pf);
    
    ///check if were still in same compartment
    coordinate = MidPointCoordinate(firstBead->coordinate, secondBead->coordinate, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    std::vector<Cylinder*> neighbors;
    ///If not initaliziation, get a neighbors list
    if(creation || extensionFront || extensionBack)  neighbors = findNearbyCylinders();
    
#ifdef CHEMISTRY
    _cCylinder = std::unique_ptr<CCylinder>(
        ChemInitializer::createCCylinder(ChemInitializerCylinderKey(), pf, _compartment, extensionFront, extensionBack, creation));
    _cCylinder->setCylinder(this);
    
    std::vector<CCylinder*> cNeighbors(neighbors.size());
    std::transform(neighbors.begin(),neighbors.end(),cNeighbors.begin(), [](Cylinder *c){return c->getCCylinder();});
    
    ///Update filament reactions, only if not initialization
    ChemInitializer::updateCCylinder(ChemInitializerCylinderKey(), _cCylinder.get(), cNeighbors);
    
#endif

#ifdef MECHANICS
    double eqLength;
    
    ///set eqLength according to cylinder size
    if(extensionFront || extensionBack) eqLength = SystemParameters::Geometry().monomerSize;
    else if(creation) eqLength = SystemParameters::Geometry().monomerSize * 2;
    else eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = std::unique_ptr<MCylinder>(new MCylinder(eqLength));
    
    _mCylinder->setCylinder(this);
    _mCylinder->setCoordinate(coordinate);
    
    ///Update neighbors list
    std::vector<MCylinder*> mNeighbors(neighbors.size());
    std::transform(neighbors.begin(),neighbors.end(),mNeighbors.begin(), [](Cylinder *c){return c->getMCylinder();});
    
    _mCylinder->updateExVolNeighborsList(mNeighbors);
#endif
    
    ///add to compartment
    _compartment->addCylinder(this);
}

Cylinder::~Cylinder() {

    ///remove from compartment
    _compartment->removeCylinder(this);
}

bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}


std::vector<Cylinder*> Cylinder::findNearbyCylinders() {
    
    std::vector<Cylinder*> cylinders;
    
    ///Find surrounding compartments (For now its conservative, change soon)
    std::vector<Compartment*> compartments;
    GController::findCompartments(coordinate, _compartment,
        SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments)
        for(auto &cyl : c->getCylinders())
            if(_pFilament != cyl->getFilament()) cylinders.push_back(cyl);

    return std::vector<Cylinder*>(cylinders.begin(), cylinders.end());
}


void Cylinder::updatePosition() {

    ///check if were still in same compartment, set new position
    coordinate = MidPointCoordinate(_pFirst->coordinate, _pSecond->coordinate, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
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
    std::vector<CCylinder*> cNeighbors(neighbors.size());
    std::transform(neighbors.begin(),neighbors.end(),cNeighbors.begin(), [](Cylinder *c){return c->getCCylinder();});
    
    ChemInitializer::updateCCylinder(ChemInitializerCylinderKey(), _cCylinder.get(), cNeighbors);
#endif
    
    ///update exvol neighbors list
#ifdef MECHANICS
    _mCylinder->setCoordinate(coordinate);
    
    std::vector<MCylinder*> mNeighbors(neighbors.size());
    std::transform(neighbors.begin(),neighbors.end(),mNeighbors.begin(), [](Cylinder *c){return c->getMCylinder();});
    
    _mCylinder->updateExVolNeighborsList(mNeighbors);
#endif
    
}


