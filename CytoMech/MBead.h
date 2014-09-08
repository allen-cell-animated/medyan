//
//  MBead.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MBead__
#define __CytoMech__MBead__

#include <iostream>
#include <vector>
#include <list>
#include "Mcommon.h"
#include "MComponent.h"
#include "BoundaryElement.h"
#include "MathFunctions.h"

using namespace mathfunc;

class Bead : public MComponent
{
public:
    
    ///The basic and simplest mecanical class. Contains information about a bead and some mecanical constant needed to calculate interactions. Constants are local and related to a bond between i(curent) and i-1 bead.
    Bead (std::vector<double> v);
    
    Bead(): coordinate (3, 0), _parent(NULL), force(3, 0), forceAux(3, 0) {}
//    virtual ~Bead();
    
	std::vector<double> coordinate; //Coordinates of the bead
	std::vector<double> force;      // Forces based on curent coordinates. Forces should alwais corespond to current coordinates.
    std::vector<double> forceAux; //An Aux field neede during CG minimization.

    //Aux functios:
    double CalcForceSquare() {return force[0]*force[0] + force[1]*force[1] + force[2]*force[2];  } // Aux method for CG minimization;
    double CalcForceSquare( int i) {return forceAux[0]*forceAux[0] + forceAux[1]*forceAux[1] + forceAux[2]*forceAux[2];  }
    double CalcDotForceProduct() { return force[0]*forceAux[0] + force[1]*forceAux[1] + force[2]*forceAux[2];}
    
    //Function which returns a pointer to a filament(parent "Composite pattern") to which this bead belongs:
    MComposite* GetParent(){
        return _parent;
    }
    
    void SetParent(MComposite* pc) {_parent = pc;}
    
    ///Set and getcompartment
    Compartment* getCompartment() {return _compartment;}
    
    ///Set the current compartment that this boundary element is in
    void setCompartment() {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBead(this);
        
        ///Add to new compartment
        _compartment = GController::getCompartment(coordinate); _compartment->addBead(this);
    }
    
    ///Alternate set compartment when compartment is known
    void setCompartment(Compartment* c) {
        
        ///remove from old compartment
        if(_compartment != nullptr) _compartment->removeBead(this);
        
        ///add to new compartment
        _compartment = c; _compartment->addBead(this);
    }
    
    ///update the boundary elements that interact with this bead
    void updateBoundaryElements();

    
private:
    
    MComposite *_parent; ///< A pointer to a parent class (Filament, linker, motor)
    Compartment* _compartment; ///< ptr to the compartment that this bead is in
    std::vector<BoundaryElement*> _boundaryElements; ///<list of currently interacting boundary elements
};


#endif /* defined(__CytoMech__MBead__) */
