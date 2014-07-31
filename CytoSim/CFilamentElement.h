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

namespace chem {
    
    class Compartment;
    class Species;
    
    /// CFilamentElement class represents a container template for all species that could be contained in a
    /// particular filament element at a given position.
    /*!
     *  CFilamentElement provides a container to hold all species that are possibly held at a given
     *  filament position. The species are held in an standard vector. Functions to lookup species
     *  as well as a filament element checker are provided.
     */
    class CFilamentElement {
        
    protected:
        Compartment* _compartment; ///< compartment that this filament element is in
        
    public:
        ///Constructor does nothing
        CFilamentElement(Compartment* c) : _compartment(c) {};
        
        ///Default destructor, removes species from compartment
        virtual ~CFilamentElement () {};
        
        ///Print a species in this filament element
        virtual void print() = 0;
        
        ///Check if this filament element is valid. Involves checking copy numbers
        virtual bool checkSpecies(int sum) = 0;
        
        ///Get species by name
        virtual Species* getSpeciesByName(std::string& name) = 0;
        
        ///move all species (and associated reactions) in this element to another compartment
        virtual void moveToCompartment(Compartment *c) = 0;
        
    };
    
    ///CMonomer class is an implementation of the abstract class CFilamentElement for a CMonomer in filament
    class CMonomer : public CFilamentElement
    {
        
    public:
        ///Constructor takes any number of species
        CMonomer(Compartment* c) : CFilamentElement(c){}
        
        ///Default destructor, does nothing
        virtual ~CMonomer () {}
        
        ///Get active filament species from this CMonomer
        virtual Species* getActiveFilamentSpecies() = 0;
        ///Get active end species from this CMonomer
        virtual Species* getActiveEndSpecies() = 0;
        
    };
    
    ///CBound class is an implementation of the abstract class CFilamentElement for a CBound on a filament
    class CBound : public CFilamentElement
    {
        
    public:
        ///Constructor takes any number of species
        CBound(Compartment* c) : CFilamentElement(c){}
        
        ///Default destructor, does nothing
        virtual ~CBound () {}
    };
}

#endif /* defined(__CytoSim__CFilamentElement__) */
