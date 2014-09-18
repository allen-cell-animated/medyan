//
//  CFilamentElement.h
//  CytoSim
//
//  Created by James Komianos on 7/17/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CFilamentElement__
#define __CytoSim__CFilamentElement__

#include <iostream>

#include "common.h"

class SpeciesFilament;
class SpeciesBound;

/// CFilamentElement class represents a container template for all species that could be contained in a
/// particular filament element at a given position.
/*!
 *  CFilamentElement provides a container to hold all species that are possibly held at a given
 *  filament position. The species are held in an standard vector. Functions to lookup species
 *  as well as a filament element checker are provided.
 */
class CFilamentElement {
    
public:
    ///Constructor does nothing
    CFilamentElement() {};
    
    ///Default destructor, removes species from compartment
    virtual ~CFilamentElement () {};
    
    ///Print a species in this filament element
    virtual void print() = 0;
    
    ///Check if this filament element is valid. Involves checking copy numbers
    virtual bool checkSpecies(int sum) = 0;
    
    ///Get species by name
    virtual Species* getSpeciesByName(std::string& name) = 0;
    
};

///CMonomer class is an implementation of the abstract class CFilamentElement for a CMonomer in filament
class CMonomer : public CFilamentElement
{
protected:
    std::vector<SpeciesFilament*> _species; ///Vector of filament species that this monomer contains
    
public:
    ///Constructor takes any number of species
    CMonomer() {}
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new compartment. The
    /// original species will remain intact, while new identical species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c) {
        
        for(auto &s: rhs._species) {
            SpeciesFilament* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            _species.push_back(sNew);
        }
    }
    
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    ///Move constructor, simply copies species vector
    CMonomer(CMonomer &&rhs) noexcept : _species(rhs._species) {};
    
    ///Move assigment operator, same as move constructor
    CMonomer& operator=(CMonomer&& rhs)  {
        _species = rhs._species;
        return *this;
    }
    
    ///Default destructor, does nothing
    virtual ~CMonomer () {}
    
    ///Clone, calls copy constructor
    virtual CMonomer* clone(Compartment* c) = 0;
    
    ///Get the vector of species
    std::vector<SpeciesFilament*>& species() {return _species;}
    
    ///Get species by name
    virtual Species* getSpeciesByName(std::string& name) {
        
        for (auto &s : _species)
            if(name.find(s->getName()) != std::string::npos) return s;
        return nullptr;
    }

    
};

///CBound class is an implementation of the abstract class CFilamentElement for a CBound on a filament
class CBound : public CFilamentElement
{
protected:
    std::vector<SpeciesBound*> _species; ///Vector of filament species that this bound contains
    
public:
    ///Constructor takes any number of species
    CBound() {}
    
    /// Copy constructor
    /// This constructor will create a new CBound, identical to the copied, in a new compartment. The
    /// original species will remain intact, while new identical species will be initialized.
    CBound(const CBound& rhs, Compartment* c) {
        
        for(auto &s: rhs._species) {
            SpeciesBound* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            _species.push_back(sNew);
        }
    }
    
    /// Assignment is not allowed
    CBound& operator=(CBound &rhs) = delete;
    
    ///Move constructor, simply copies species vector
    CBound(CBound &&rhs) noexcept : _species(rhs._species) {};
    
    ///Move assigment operator, same as move constructor
    CBound& operator=(CBound&& rhs)  {
        _species = rhs._species;
        return *this;
    }
    ///Default destructor, does nothing
    virtual ~CBound () {}
    
    ///Clone, calls copy constructor
    virtual CBound* clone(Compartment* c) = 0;
    
    ///Get the vector of species
    std::vector<SpeciesBound*>& species() {return _species;}
    
    ///Get species by name
    virtual Species* getSpeciesByName(std::string& name) {
        
        for (auto &s : _species)
            if(name.find(s->getName()) != std::string::npos) return s;
        return nullptr;
    }
    
};


#endif /* defined(__CytoSim__CFilamentElement__) */
