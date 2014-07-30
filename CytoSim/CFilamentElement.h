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
#include "Compartment.h"
#include "Species.h"

namespace chem {
    
    /// CFilamentElement class represents a container template for all species that could be contained in a
    /// particular filament element at a given position.
    /*!
     *  CFilamentElement provides a container to hold all species that are possibly held at a given
     *  filament position. The species are held in an standard array. Functions to lookup species
     *  as well as a filament element checker are provided.
     */
    class CFilamentElement {
        
    protected:
        Compartment* _compartment; ///< compartment that this filament element is in
        
    public:
        ///Constructor does nothing
        CFilamentElement(Compartment* c) : _compartment(c) {};
        
        ///Default destructor, removes species from compartment
        ~CFilamentElement ()
        {
            //            for(auto& s : _species) {
            //                //_compartment->removeInternalReactions(s);
            //                //_compartment->removeSpecies(s);
            //            }
        };
        
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
        ~CMonomer () {}
        
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
        ~CBound () {}
    };
    
    ///Basic monomer consisting of actin, formin, capping, and a virtual front/back
    class CMonomerBasic : public CMonomer {
        
    private:
        SpeciesFilament* _filament_species[1]; ///<array of filament species (just actin for now)
        SpeciesFilament* _end_species[4]; ///<array of "end" species (front, back, formin, capping)
        
    public:
        ///Constructor, initializes species container
        CMonomerBasic(std::vector<SpeciesFilament*> species, Compartment* c);
        
        ///Destructor, removes all species and associated reactions from compartment
        ~CMonomerBasic();
        
        ///Get a species
        virtual SpeciesFilament* getActin() {return _filament_species[0];}
        
        virtual SpeciesFilament* getFront() {return _end_species[0];}
        
        virtual SpeciesFilament* getBack() {return _end_species[1];}
        
        virtual SpeciesFilament* getFormin() {return _end_species[2];}
        
        virtual SpeciesFilament* getCapping() {return _end_species[3];}
        
        ///Look up species by name
        virtual Species* getSpeciesByName(std::string& name);
        
        
        ///Check if this monomer is valid
        virtual bool checkSpecies(int sum)
        {
            return true;
//            int currentSum = 0;
//            for(auto &s : _species)
//                currentSum += s->getN();
//            return currentSum = sum;
        }

        ///Print a species in this filament element
        virtual void print();
        
        ///move all species (and associated reactions) in this element to another compartment
        virtual void moveToCompartment(Compartment *c) { 
        }
        
        ///Find active filament species
        ///@note return null if none
        virtual Species* getActiveFilamentSpecies();
        
        ///Find active end species
        ///@note return null if none
        virtual Species* getActiveEndSpecies();
    };
    
    ///Basic bound consisting of myosin, myosin-actin, and a virtual empty species
    class CBoundBasic : public CBound {
        
    private:
        SpeciesBound* _species[3]; ///<array of contained species
        
    public:
        CBoundBasic(std::vector<SpeciesBound*> species, Compartment* c) : CBound(c);
        
        ///Destructor, removes all species and associated reactions from compartment
        ~CBoundBasic();
        
        ///Look up a species given a name
        virtual SpeciesBound* getEmpty() {return _species[0];}
        
        virtual SpeciesBound* getMyosin() {return _species[1];}
        
        virtual SpeciesBound* getMyosinActin() {return _species[2];}
        
        ///Check if this monomer is valid
        virtual bool checkSpecies(int sum);
        
        ///Look up species by name
        virtual Species* getSpeciesByName(std::string& name);
        
        ///Print a species in this filament element
        virtual void print();

        virtual void moveToCompartment(Compartment* c) {}
    };
    
}

#endif /* defined(__CytoSim__CFilamentElement__) */
