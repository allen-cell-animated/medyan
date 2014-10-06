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
#include "Compartment.h"
#include "MathFunctions.h"
#include "GController.h"


using namespace mathfunc;

///The Bead class represents a single coordinate and mechanical constants needed.
/*!
 * The basic and simplest mechanical class. Contains information about a bead and some mecanical
 * constant needed to calculate interactions. Constants are local and related to a bond between i(curent)
 * and i-1 bead.
 */


class Bead : public Component {
public:

    std::vector<double> coordinate; ///< Coordinates of the bead
	std::vector<double> force;
    ///< Forces based on curent coordinates. Forces should always correspond to current coordinates.
    std::vector<double> forceAux; ///< An Aux field neede during CG minimization.
    
    ///Main constructor
    Bead (std::vector<double> v, int ID);
    ///Default constructor
    Bead(int ID): coordinate (3, 0), force(3, 0), forceAux(3, 0), _ID(ID) {}

    ///Aux functions
    
    // Aux method for CG minimization
    double CalcForceSquare() {return force[0]*force[0] + force[1]*force[1] + force[2]*force[2]; }
    double CalcForceSquare(int i) {return forceAux[0]*forceAux[0] + forceAux[1]*forceAux[1] + forceAux[2]*forceAux[2];}
    double CalcDotForceProduct() { return force[0]*forceAux[0] + force[1]*forceAux[1] + force[2]*forceAux[2];}
    
    ///Set and get compartment
    Compartment* getCompartment() {return _compartment;}
    
    ///Set the current compartment that this boundary element is in
    void setCompartment() {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBead(this);
        
        ///Add to new compartment
        try {_compartment =  GController::getCompartment(coordinate);}
        catch (std::exception& e) {std::cout << e.what(); exit(EXIT_FAILURE);}
        _compartment->addBead(this);
    }
    
    ///Alternate set compartment when compartment is known
    void setCompartment(Compartment* c) {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBead(this);
        ///add to new compartment
        _compartment = c; _compartment->addBead(this);
    }
    
    ///update the position of this bead (and interaction boundary elements)
    void updatePosition();
    
    ///getters for bead data
    int getID() {return _ID;}
    float getBirthTime() {return _birthTime;}

    
private:
    Compartment* _compartment = nullptr; ///< ptr to the compartment that this bead is in
    std::vector<BoundaryElement*> _boundaryElements; ///<list of currently interacting boundary elements
    
    int _ID; ///Unique ID of the bead in this filament (relative to the filament)
    float _birthTime; ///Time of birth of bead;
};


#endif /* defined(__CytoMech__Bead__) */
