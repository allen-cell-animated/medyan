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
    ///ALL SPECIES VECTORS
    std::vector<SpeciesFilament*> _speciesFilament; ///<Filament species
    std::vector<SpeciesBound*> _speciesBound;       ///<Bound species
    std::vector<SpeciesPlusEnd*> _speciesPlusEnd;   ///<PlusEnd species
    std::vector<SpeciesMinusEnd*> _speciesMinusEnd; ///<MinusEnd species
public:
    ///Constructor does nothing
    CMonomer() {};
    
    ///Default destructor
    virtual ~CMonomer () {};
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new compartment. The
    /// original species will remain intact, while new identical species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c);
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    ///Move constructor, simply copies species vectors
    CMonomer(CMonomer &&rhs) noexcept : _speciesFilament(rhs._speciesFilament),
                                        _speciesBound(rhs._speciesBound),
                                        _speciesPlusEnd(rhs._speciesPlusEnd),
                                        _speciesMinusEnd(rhs._speciesMinusEnd) {};
    
    ///Move assigment operator, same as move constructor
    CMonomer& operator=(CMonomer&& rhs)  {
        _speciesFilament = rhs._speciesFilament;
        _speciesBound = rhs._speciesBound;
        
        _speciesPlusEnd = rhs._speciesPlusEnd;
        _speciesMinusEnd = rhs._speciesMinusEnd;
        return *this;
    }
    
    virtual CMonomer* clone(Compartment* c) {
        return new CMonomer(*this, c);
    }
    
    ///Add a species filament
    void addSpeciesFilament(SpeciesFilament* s) { _speciesFilament.push_back(s); }
    ///Add a species bound
    void addSpeciesBound(SpeciesBound* s) { _speciesBound.push_back(s); }
    ///Add a species minus end
    void addSpeciesPlusEnd(SpeciesPlusEnd* s) { _speciesPlusEnd.push_back(s); }
    ///Add a species plus end
    void addSpeciesMinusEnd(SpeciesMinusEnd* s) { _speciesMinusEnd.push_back(s); }
    
    ///Print the species in this filament element
    void print();

    ///Return all species vectors
    const std::vector<SpeciesFilament*>& speciesFilamentVector() {return _speciesFilament;}
    const std::vector<SpeciesBound*>& speciesBoundVector() {return _speciesBound;}
    const std::vector<SpeciesPlusEnd*>& speciesPlusEndVector() {return _speciesPlusEnd;}
    const std::vector<SpeciesMinusEnd*>& speciesMinusEndVector() {return _speciesMinusEnd;}
    
    ///Get species at a specific index
    ///@note no check on this index. The index value of a species is stored in the chemical initializer
    ///when all reactions are initialized from the chemical input file
    SpeciesFilament* speciesFilament(int index) {return _speciesFilament[index];}
    SpeciesBound* speciesBound(int index) {return _speciesBound[index];}
    SpeciesPlusEnd* speciesPlusEnd(int index) {return _speciesPlusEnd[index];}
    SpeciesMinusEnd* speciesMinusEnd(int index) {return _speciesMinusEnd[index];}
    
    ///Check if this filament element is valid. Involves checking copy numbers
    virtual bool checkSpecies(int sum) {return true;}
    
};




#endif /* defined(__CytoSim__CFilamentElement__) */
