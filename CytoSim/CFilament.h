//
//  CFilament.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CFilament__
#define __CytoSim__CFilament__

#include <iostream>
#include "CompartmentContainer.h"


namespace chem {
    
    /// CFilamentElement class represents a container template for all species that could be contained in a
    /// particular filament element at a given position.
    /*!
     *  CFilamentElement provides a container to hold all species that are possibly held at a given
     *  filament position. The species are held in an standard array. Functions to lookup species
     *  as well as a filament element checker are provided.
     */
    template <class SpeciesType>
    class CFilamentElement {
        
    protected:
        std::vector<SpeciesType*> _species; ///< array of species in this element
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
        
        ///Add a species to this CFilamentElement
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
    
    ///Monomer class is an implementation of the abstract class CFilamentElement for a monomer in filament
    class Monomer : public CFilamentElement<SpeciesFilament>
    {
        
    public:
        ///Constructor takes any number of species
        Monomer(std::vector<SpeciesFilament*> species, Compartment* c) :
            CFilamentElement<SpeciesFilament>(c)
        {
            for(auto &s: species){CFilamentElement<SpeciesFilament>::_species.push_back(s);}
        }
        
        ///Default destructor, does nothing
        ~Monomer () {}
        
        ///Check if this monomer is valid.
        virtual bool checkSpecies(int sum)
        {
            int currentSum = 0;
            for(auto &s : CFilamentElement<SpeciesFilament>::_species) {currentSum += s->getN();}
            return currentSum = sum;
        }
        
    };
    
    ///Bound class is an implementation of the abstract class CFilamentElement for a monomer in filament
    class Bound : public CFilamentElement<SpeciesBound>
    {
        
    public:
        ///Constructor takes any number of species
        Bound(std::vector<SpeciesBound*> species, Compartment* c) :
            CFilamentElement<SpeciesBound>(c)
        {
            for(auto &s: species) {CFilamentElement<SpeciesBound>::_species.push_back(s);}
        }
        
        ///Default destructor, does nothing
        ~Bound () {}
        
        ///Check if this monomer is valid.
        virtual bool checkSpecies(int sum)
        {
            int currentSum = 0;
            for(auto &s : CFilamentElement<SpeciesBound>::_species) {currentSum += s->getN();}
            return currentSum = sum;
        }
    };
    
    
    /// CSubFilament class holds all Monomers and Bounds
    /*! 
     *  The CSubFilamentClass is an template class that has lists of the monomers and bounds that it contains.
     *  it has functionality to print the current composition.
     *  Accessing a particular species in the CSubFilament is possible as well.
     */
    class CSubFilament : public Component {
        
    protected:
        std::vector<std::unique_ptr<Monomer>> _monomers; ///< list of monomers in this sub filament
        std::vector<std::unique_ptr<Bound>> _bounds; ///< list of bound species in this sub filament
        Compartment* _compartment; ///< compartment this CSubFilament is in
        short _full_length = 0; ///< length of this CSubFilament
        short _length = 0; 
        
    public:
        ///Default constructor, sets compartment
        CSubFilament(Compartment* c) : _compartment(c) {}
        
        ///Default destructor, explicitly removes monomers and bounds
        ~CSubFilament()
        {
            _monomers.clear();
            _bounds.clear();
        }
        
        ///get filament compartment
        Compartment* compartment() {return _compartment;}
        
        ///Add a monomer to this CSubFilament
        virtual void addMonomer(Monomer* monomer) {
            _monomers.emplace_back(std::unique_ptr<Monomer>(monomer));
            _full_length++;
        }
        
        ///Add a bound to this CSubFilament
        virtual void addBound(Bound* bound) {_bounds.emplace_back(std::unique_ptr<Bound>(bound));}
        
        ///Get monomer at an index
        ///@note no check on index
        virtual Monomer* monomer(int index) {return _monomers[index].get();}
        
        ///Get bound at an index
        ///@note no check on index
        virtual Bound* bound(int index) {return _bounds[index].get();}
        
        ///Get back or front monomer/bound
        ///@note monomer/bound list must not be empty
        virtual Monomer* backMonomer() {return _monomers[_full_length - 1].get();}
        virtual Bound* backBound() {return _bounds[_full_length - 1].get();}
        virtual Monomer* frontMonomer() {return _monomers[0].get();}
        virtual Bound* frontBound() {return _bounds[0].get();}
        
        ///Get species at specified index (monomer)
        virtual SpeciesFilament* getMonomerSpecies(int index, std::string name) {
            
            return monomer(index)->species(name);
        }
        
        ///Get species at specified index (bound)
        virtual SpeciesBound* getBoundSpecies(int index, std::string name) {
            
            return bound(index)->species(name);
        }
        
        ///Get the current length
        virtual short length() {return _length;}
        ///Increase length
        virtual void increaseLength() {if(_length != _full_length) _length++;}
        ///Decrease length
        virtual void decreaseLength() {if(_length != 0) _length--;}
        ///see if the subfilament is at maxlength
        virtual bool atMaxLength() {return _length == _full_length;}
        ///set the length of this subfilament
        virtual void setLength(int length) {_length = length;}
        
        
        ///Print CSubFilament
        virtual void printCSubFilament()
        {
            std::cout << "Composition of CSubFilament: " << std::endl;
            for (auto &m : _monomers){
                m->print();
                std::cout << ":";
            }
            std::cout << std::endl << "Bounds of CSubFilament: " <<std::endl;
            for (auto &b : _bounds) {
                b->print();
                std::cout << ":";
            }
        }
    };
    
    
    /// CFilament class holds sub CFilaments
    /*! The CFilament class is used to hold CSubFilaments (children of CFilament).
     */
    class CFilament : public Composite{
    
    private:
        int _length = 0; ///Length of filament
        
    public:
        ///Default constructor, does nothing
        CFilament() {};
        
        ///Default destructor, removes all Csubfilaments implicitly
        ~CFilament() {}
        
        ///Add a Csubfilament
        virtual void addCSubFilament(CSubFilament* s) {
            addChild(std::unique_ptr<Component>(s));
        }
        
        ///Get front subfilament
        ///@note -  no check on the number of children
        virtual CSubFilament* getFrontCSubFilament()
        {
            CSubFilament* front = nullptr;
            int childIndex = numberOfChildren() - 1;
            
            while(childIndex >= 0) {
                front = static_cast<CSubFilament*>(children(childIndex));
                
                if (front->length() != 0) return front;
                else {
                    if((childIndex - 1) >= 0) {
                        CSubFilament* second =
                            static_cast<CSubFilament*>(children(childIndex - 1));
                        if (second->atMaxLength())
                            return front;
                    }
                }
                childIndex--;
            }
            return front;
        }
        
        ///number of subfilaments in this filament
        virtual int numCSubFilaments() {return numberOfChildren();}
        
        ///Increase length of filament
        virtual void increaseLength() {
            getFrontCSubFilament()->increaseLength();
            _length++;
        }
        ///Decrease length of filament
        virtual void decreaseLength() {
            getFrontCSubFilament()->decreaseLength();
            _length--;
        }
        
        ///Set the length of this filament
        ///@note should only be used when the filament has EXACTLY ONE sub filament
        virtual void setLength(int length) {
            getFrontCSubFilament()->setLength(length);
            _length = length;
        }
        
        ///get length of filament
        virtual int length() {return _length;}
        
        ///Print entire filament
        virtual void printCFilament() {
            
            std::cout<<std::endl;
            std::cout << "Length = " << _length <<std::endl;
            int index = 0;
            for (auto &c : children()) {
                std::cout << "CSubFilament " << index++ << ":" <<std::endl;
                static_cast<CSubFilament*>(c.get())->printCSubFilament();
                std::cout<< std::endl<<std::endl;
            }
        }
        
    };
    

    /// Membrane class updates polymerization rates based on the Brownian-Ratchet model
    
    /*! The Membrane class will update all polymerization reaction rates based on
     *  the given configuration of filament ends and a membrane. The calculation is:
     *
     *
     *  1) The effective polymerization rate of a filament is equal to the "bare" rate,
     *     denoted as k0, times the probability of the gap between the filament and the
     *     membrane opening. The relation between the loading force and the rate is:
     *
     *                      k_n + = k+ * exp(-f_n * monomer size / k * T)
     *
     *     where f_n is the current force on the nth filament.
     *
     *  2) To compute the forces f_n, we must calculate the probability of finding the
     *     nth filament below the membrane (which we denote as a height of h = max_n(h_n))
     *
     *             p_n ~ integral from h-h_n to infinity ( exp(-z^2 / sigma^2) dz )
     *
     *     where sigma is the average membrane fluctuation amplitude, which we assume
     *     to be constant.
     *
     *  3) Then, the forces are computed as f_n = (p_n / p) * f, where the total force f
     *     is constant, and p is the normalization factor.
     *
     */
    
    class Membrane {
        
    private:
        std::vector<ReactionBase*> _poly_reactions; ///<vector of reactions to update
        
        int _h_membrane = 0; ///<height of the membrane
        
        
        
        
        
        
    };
    
    
    
    
    
    
    
    
    
    
}; //namespace chem








#endif /* defined(__CytoSim__CFilament__) */
