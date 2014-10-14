//
//  CBound.h
//  Cyto
//
//  Created by James Komianos on 10/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CBound__
#define __Cyto__CBound__

#include <iostream>
#include "Species.h"

class Compartment;

///CBound class represents a chemical object that is bound to a filament
/*!
 *  The CBound class is an abstract representation of a chemically bound object to a filament
 *  (Could be a linker, motor, branching protein, etc). Each CBound object has a pointer
 *  to the corresponding species on a filament. Different implementations of CBound class will
 *  have different functions to bind, move, etc.
 */
class CBound {
    
protected:
    SpeciesBound* _firstSpecies = nullptr; ///< corresponding species on filament
    SpeciesBound* _secondSpecies = nullptr; ///< second species the CBound could bind to on another filament
    
    Compartment* _compartment; ///< compartment this CBound is in
    
public:
    ///Constructor, just sets species
    CBound(Compartment* c) : _compartment(c) {}
    
    ///virtual destructor
    virtual ~CBound() noexcept {
        _firstSpecies->removeCBound();
        _secondSpecies->removeCBound();
    }
    
    ///setter and getter for first species
    void setFirstSpecies(SpeciesBound* species) {
        ///remove from old first species
        if(_firstSpecies != nullptr) _firstSpecies->removeCBound();
        
        _firstSpecies = species;
        _firstSpecies->setCBound(this);
    }
    SpeciesBound* getFirstSpecies() {return _firstSpecies;}
    
    ///setter and getter for second species
    void setSecondSpecies(SpeciesBound* species) {
        ///remove from old second species
        if(_secondSpecies != nullptr) _secondSpecies->removeCBound();
        
        _secondSpecies = species;
        _secondSpecies->setCBound(this);
    }
    SpeciesBound* getSecondSpecies() {return _secondSpecies;}
    
    ///getter for compartment
    Compartment* getCompartment() {return _compartment;}
};

#endif /* defined(__Cyto__CBound__) */
