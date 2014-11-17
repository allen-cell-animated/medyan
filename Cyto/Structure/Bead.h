//
//  Bead.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__Bead__
#define __CytoMech__Bead__

#include <iostream>
#include <vector>
#include <list>

#include "common.h"

#include "Component.h"
#include "Neighbor.h"
#include "Compartment.h"
#include "GController.h"

///FORWARD DECLARATIONS
class Compartment;
class BoundaryElement;

///The Bead class represents a single coordinate and mechanical constants needed.
/*!
 * The basic and simplest mechanical class. Contains information about a bead and some mecanical
 * constant needed to calculate interactions. Constants are local and related to a bond between i(curent)
 * and i-1 bead.
 */

class Bead : public Component, public Neighbor {
public:

    vector<double> coordinate; ///< Coordinates of the bead
    vector<double> coordinateAux; ///< An aux coordinate field needed during CG minimization
	vector<double> force;
    ///< Forces based on curent coordinates. Forces should always correspond to current coordinates.
    vector<double> forceAux; ///< An Aux field needed during CG minimization.
    
    ///Main constructor
    Bead (vector<double> v, int positionFilament);
    ///Default constructor
    Bead(int positionFilament) : _positionFilament(positionFilament),
                                 coordinate (3, 0), coordinateAux(3, 0),
                                 force(3, 0), forceAux(3, 0) {}
    ~Bead();
    
    ///Aux functions
    
    // Aux method for CG minimization
    double calcForceSquare() {return force[0]*force[0] + force[1]*force[1] + force[2]*force[2]; }
    double calcForceAuxSquare() {return forceAux[0]*forceAux[0] + forceAux[1]*forceAux[1] + forceAux[2]*forceAux[2];}
    double calcDotForceProduct() { return force[0]*forceAux[0] + force[1]*forceAux[1] + force[2]*forceAux[2];}
    
    ///add a boundary element to list of interacting boundary elements
    void addBoundaryElement(BoundaryElement* be) {_boundaryElements.push_back(be);}
    ///Remove a boundary element from list of interacting boundary elements
    ///@note does nothing if boundary element is not in interacting list already
    void removeBoundaryElement(BoundaryElement* be) {
        auto it = find(_boundaryElements.begin(), _boundaryElements.end(), be);
        if(it != _boundaryElements.end()) _boundaryElements.erase(it);
    }
    
    ///Set and get compartment
    Compartment* getCompartment() {return _compartment;}
    
    ///update the position of this bead (and interaction boundary elements)
    void updatePosition();
    
    ///getters for bead data
    void setPositionFilament(int positionFilament) {_positionFilament = positionFilament;}
    int getPositionFilament() {return _positionFilament;}
    
    float getBirthTime() {return _birthTime;}

    
private:
    Compartment* _compartment = nullptr; ///< ptr to the compartment that this bead is in
    vector<BoundaryElement*> _boundaryElements; ///<list of currently interacting boundary elements
    
    int _positionFilament; ///Position of bead on filament
    float _birthTime; ///Time of birth of bead;
};


#endif /* defined(__CytoMech__Bead__) */
