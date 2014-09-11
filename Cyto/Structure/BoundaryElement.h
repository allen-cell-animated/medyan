//
//  BoundaryElement.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElement__
#define __Cyto__BoundaryElement__

#include <vector>
#include <iostream>
#include "Gcontroller.h"

class Bead;

///BoundaryElement class represents an element of a boundary surface
/*!
 * The BoundaryElement class is a representation of a boundary surface element, which can interact
 * with other elements in the system, including other BoundaryElements as well as Beads in filaments.
 * Together, a collection of boundary elements make up a BoundarySurface.
 */
class BoundaryElement {
    
protected:
    
    Compartment* _compartment = nullptr; ///< Compartment that this boundary element is currently in
    std::vector<Bead*> _beads; ///< Beads that this boundary element could interact with
    std::vector<BoundaryElement*> _neighbors; ///<neighbors of this boundary element
    
    std::vector<double> _coords; ///< coordinates of this boundary element
    
public:
    ///Default constructor
    BoundaryElement(std::vector<double> coords) : _coords(coords) {
    
        ///set the compartment given the initial coordinates
        setCompartment();
    }
    
    ///Destructor
    ~BoundaryElement() {
        
        ///remove from compartment
        _compartment->removeBoundaryElement(this);
        
        ///remove from neighbors
        for (auto &be : _neighbors) be->removeNeighbor(this);
    }
    
    ///add a bead to list of interacting beads
    void addBead(Bead* b) {_beads.push_back(b);}
    ///Remove a bead from list of interacting beads
    ///@note does nothing if bead is not in interacting list already
    void removeBead(Bead* b) {
        auto it = std::find(_beads.begin(), _beads.end(), b);
        if(it != _beads.end()) _beads.erase(it);
    }
    
    ///Add a boundary element neighbor to this element
    void addNeighbor(BoundaryElement* b) {_neighbors.push_back(b);}
    ///remove a boundary element neighbor
    ///@note does nothing if boundary element is not in list already
    void removeNeighbor(BoundaryElement* b) {
        auto it = std::find(_neighbors.begin(), _neighbors.end(), b);
        if(it != _neighbors.end()) _neighbors.erase(it);
    }
    
    ///Set the current compartment that this boundary element is in
    void setCompartment() {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBoundaryElement(this);
        
        ///Add to new compartment
        _compartment = GController::getCompartment(_coords);
        _compartment->addBoundaryElement(this);
    }
    
    ///Alternate set compartment when compartment is known
    void setCompartment(Compartment* c) {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBoundaryElement(this);
        
        ///add to new compartment
        _compartment = c;
        _compartment->addBoundaryElement(this);
    }
    
    ///Get the compartment that this element is in
    Compartment* getCompartment() {return _compartment;}
    
    
    ///return coordinates of boundary element
    const std::vector<double>& coords() {return _coords;}
    ///set coordinates of boundary element
    void setCoords(std::vector<double>& newCoords) {_coords = newCoords;}
 
};



#endif /* defined(__Cyto__BoundaryElement__) */
