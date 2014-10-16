//
//  CMonomer.h
//  CytoSim
//
//  Created by James Komianos on 7/17/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CMonomer__
#define __CytoSim__CMonomer__

#include <iostream>

#include "common.h"
#include "Species.h"
#include "Compartment.h"

/// CMonomer class represents a container template for all species that could be contained in a
/// particular filament element at a given position.
/*!
 *  CMonomer provides a container to hold all species that are possibly held at a given
 *  filament position. The species are held in an standard vector. Functions to lookup species
 *  as well as a filament element checker are provided.
 */
class CMonomer {
    
    ///Species arrays
    SpeciesFilament* _speciesFilament[NUMSPECIESFILAMENT];
    SpeciesPlusEnd*  _speciesPlusEnd[NUMSPECIESPLUSEND];
    SpeciesMinusEnd* _speciesMinusEnd[NUMSPECIESMINUSEND];
    
    SpeciesBound*  _speciesBound[NUMSPECIESBOUND];
    SpeciesLinker* _speciesLinker[NUMSPECIESLINKER];
    SpeciesMotor*  _speciesMotor[NUMSPECIESMOTOR];
    
public:
    ///Constructor does nothing but reset arrays
    CMonomer();
    
    ///Default destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~CMonomer () noexcept {};
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new compartment. The
    /// original species will remain intact, while new identical species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c);
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    virtual CMonomer* clone(Compartment* c) {
        return new CMonomer(*this, c);
    }
    ///Add a species
    void addSpeciesFilament(SpeciesFilament* s);
    void addSpeciesPlusEnd(SpeciesPlusEnd* s);
    void addSpeciesMinusEnd(SpeciesMinusEnd* s);
    void addSpeciesBound(SpeciesBound* s);
    void addSpeciesLinker(SpeciesLinker* s);
    void addSpeciesMotor(SpeciesMotor* s);
    
    ///Print the species in this filament element
    void print();

    ///Get species at a specific index
    ///@note no check on this index. The index value of a species is stored in the chemical initializer
    ///when all reactions are initialized from the chemical input file
    SpeciesFilament* speciesFilament(int index) {return _speciesFilament[index];}
    SpeciesPlusEnd* speciesPlusEnd(int index) {return _speciesPlusEnd[index];}
    SpeciesMinusEnd* speciesMinusEnd(int index) {return _speciesMinusEnd[index];}
    
    SpeciesBound* speciesBound(int index) {return _speciesBound[index];}
    SpeciesLinker* speciesLinker(int index) {return _speciesLinker[index];}
    SpeciesMotor* speciesMotor(int index) {return _speciesMotor[index];}
    
    ///Check if this filament element is valid. Involves checking copy numbers
    virtual bool checkSpecies(int sum) {return true;}
    
};




#endif /* defined(__CytoSim__CMonomer__) */
