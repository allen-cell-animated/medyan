//
//  Filament.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__Filament__
#define __CytoSim__Filament__

#include <iostream>
#include "CompartmentContainer.h"


namespace chem {
    
    /// FilamentElement class represents a container template for all species that could be contained in a
    /// particular filament element at a given position.
    /*!
     *  FilamentElement provides a container to hold all species that are possibly held at a given
     *  filament position. The species are held in an standard array. Functions to lookup species
     *  as well as a filament element checker are provided.
     */
    template <class SpeciesType>
    class FilamentElement {
        
    protected:
        std::vector<SpeciesType*> _species; ///< array of species in this element
        Compartment* _compartment; ///< compartment that this filament element is in
        
    public:
        ///Constructor does nothing
        FilamentElement(Compartment* c) : _compartment(c) {};
        
        ///Default destructor, removes species from compartment
        ~FilamentElement ()
        {
            for(auto& s : _species) {
                _compartment->removeSpecies(s);
                _compartment->removeInternalReactions(s);
                s->~Species();
                delete s;
            }
        };
        
        ///Add a species to this filamentElement
        virtual void addSpecies(SpeciesType* s) {
            _species.push_back(s);
        }
        
        
        ///Look up a species given a name
        virtual SpeciesType* species(std::string name)
        {
            for (SpeciesType* &s : _species)
                if(s->getName() == name) return s;
            return nullptr;
        }
        
        ///Print a species in this filament element
        virtual void print()
        {
            for (SpeciesType* &s : _species)
                if(s->getN() == 1) std::cout << s->getName().at(0);
            
        }
        
        ///Check if this filament element is valid. Involves checking copy numbers
        virtual bool checkSpecies(int sum) = 0;
    };
    
    ///Monomer class is an implementation of the abstract class FilamentElement for a monomer in filament
    class Monomer : public FilamentElement<SpeciesFilament>
    {
        
    public:
        ///Constructor takes any number of species
        Monomer(std::vector<SpeciesFilament*> species, Compartment* c) :
            FilamentElement<SpeciesFilament>(c)
        {
            for(auto &s: species){FilamentElement<SpeciesFilament>::_species.push_back(s);}
        }
        
        ///Default destructor, does nothing
        ~Monomer () {}
        
        ///Check if this monomer is valid.
        virtual bool checkSpecies(int sum)
        {
            int currentSum = 0;
            for(auto &s : FilamentElement<SpeciesFilament>::_species) {currentSum += s->getN();}
            return currentSum = sum;
        }
        
    };
    
    ///Bound class is an implementation of the abstract class FilamentElement for a monomer in filament
    class Bound : public FilamentElement<SpeciesBound>
    {
        
    public:
        ///Constructor takes any number of species
        Bound(std::vector<SpeciesBound*> species, Compartment* c) :
            FilamentElement<SpeciesBound>(c)
        {
            for(auto &s: species) {FilamentElement<SpeciesBound>::_species.push_back(s);}
        }
        
        ///Default destructor, does nothing
        ~Bound () {}
        
        ///Check if this monomer is valid.
        virtual bool checkSpecies(int sum)
        {
            int currentSum = 0;
            for(auto &s : FilamentElement<SpeciesBound>::_species) {currentSum += s->getN();}
            return currentSum = sum;
        }
    };
    
    
    /// SubFilament class holds all Monomers and Bounds
    /*! 
     *  The SubFilamentClass is an template class that has lists of the monomers and bounds that it contains.
     *  it has functionality to print the current composition.
     *  Accessing a particular species in the subfilament is possible as well.
     */
    class SubFilament : public Component {
        
    protected:
        std::vector<std::unique_ptr<Monomer>> _monomers; ///< list of monomers in this sub filament
        std::vector<std::unique_ptr<Bound>> _bounds; ///< list of bound species in this sub filament
        Compartment* _compartment; ///< compartment this subfilament is in
        short _length = 0; ///< length of this subfilament
        
    public:
        ///Default constructor, sets compartment
        SubFilament(Compartment* c) : _compartment(c) {}
        
        ///Default destructor, implicitly removes monomers and bounds
        ~SubFilament()
        {
            _monomers.clear();
            _bounds.clear();
        }
        
        ///get filament compartment
        Compartment* compartment() {return _compartment;}
        
        ///Add a monomer to this subfilament
        virtual void addMonomer(Monomer* monomer) {
            _monomers.emplace_back(std::unique_ptr<Monomer>(monomer));
            _length++;
        }
        
        ///Add a bound to this subfilament
        virtual void addBound(Bound* bound) {_bounds.emplace_back(std::unique_ptr<Bound>(bound));}
        
        ///Get monomer at an index
        ///@note no check on index
        virtual Monomer* monomer(int index) {return _monomers[index].get();}
        
        ///Get bound at an index
        ///@note no check on index
        virtual Bound* bound(int index) {return _bounds[index].get();}
        
        ///Get end monomer
        ///@note monomer list must not be empty
        virtual Monomer* backMonomer() {return _monomers[_length - 1].get();}
        
        ///Get end monomer
        ///@note bounds list must not be empty
        virtual Bound* backBound() {return _bounds[_length - 1].get();}
        
        ///Get front monomer
        ///@note monomer list must not be empty
        virtual Monomer* frontMonomer() {return _monomers[0].get();}
        
        ///Get end monomer
        ///@note bounds list must not be empty
        virtual Bound* frontBound() {return _bounds[0].get();}
        
        ///Get species at specified index (monomer)
        virtual SpeciesFilament* getMonomerSpecies(int index, std::string name) {
            
            return monomer(index)->species(name);
        }
        
        ///Get species at specified index (bound)
        virtual SpeciesBound* getBoundSpecies(int index, std::string name) {
            
            return bound(index)->species(name);
        }
        
        ///Print subfilament
        virtual void printSubFilament()
        {
            std::cout << "Composition of SubFilament: " << std::endl;
            for (auto &m : _monomers){
                m->print();
                std::cout << ":";
            }
            std::cout << std::endl << "Bounds of SubFilament: " <<std::endl;
            for (auto &b : _bounds) {
                b->print();
                std::cout << ":";
            }
        }
    };
    
    
    /// Filament class holds sub filaments
    /*! The Filament class is used to hold subfilaments (children of filament).
     */
    class Filament : public Composite{
        
    public:
        ///Default constructor, does nothing
        Filament() {};
        
        ///Default destructor, removes all subfilaments implicitly
        ~Filament() {}
        
        ///Add a subfilament
        virtual void addSubFilament(SubFilament* s) {
            addChild(std::unique_ptr<Component>(s));
        }
        
        ///Get front subfilament
        ///@note -  no check on the number of children
        virtual SubFilament* getFrontSubFilament()
        {
            return static_cast<SubFilament*>(children(numberOfChildren() - 1));
            
        }
        
        ///number of subfilaments in this filament
        virtual int numSubFilaments() {return numberOfChildren();}
        
        ///Print entire filament
        virtual void printFilament() {
            
            std::cout<<std::endl;
            int index = 0;
            for (auto &c : children()) {
                std::cout << "SubFilament " << index++ << ":" <<std::endl;
                static_cast<SubFilament*>(c.get())->printSubFilament();
                std::cout<< std::endl<<std::endl;
            }
        }
        
    };
    

    
}; //namespace chem








#endif /* defined(__CytoSim__Filament__) */
