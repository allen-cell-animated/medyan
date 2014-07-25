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
#include "Mcommon.h"

class Filament;


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
    
    
    /// CSubFilament class holds all Monomers and Bounds
    /*! 
     *  The CSubFilament Class is an template class that has lists of the monomers and bounds that it contains.
     *  it has functionality to print the current composition.
     *  Accessing a particular species in the CSubFilament is possible as well.
     */
    class CSubFilament : public Component {
        
    protected:
        std::vector<std::unique_ptr<CMonomer>> _monomers; ///< list of monomers in this sub filament
        std::vector<std::unique_ptr<CBound>> _bounds; ///< list of bound species in this sub filament
        std::vector<ReactionBase*> _reactions;///< list of reactions associated with this subfilament
        Compartment* _compartment; ///< compartment this CSubFilament is in
        short _length = 0;
        const short _max_length = L / monomer_size;
        
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
        virtual void addCMonomer(CMonomer* monomer) {
            _monomers.emplace_back(std::unique_ptr<CMonomer>(monomer));
        }
        
        ///Add a bound to this CSubFilament
        virtual void addCBound(CBound* bound) {_bounds.emplace_back(std::unique_ptr<CBound>(bound));}
        
        ///Get monomer at an index
        ///@note no check on index
        virtual CMonomer* getCMonomer(int index) {return _monomers[index].get();}
        
        ///Get bound at an index
        ///@note no check on index
        virtual CBound* getCBound(int index) {return _bounds[index].get();}
        
        ///Get back or front monomer/bound
        ///@note monomer/bound list must not be empty
        virtual CMonomer* backCMonomer() {return _monomers[0].get();}
        
        virtual CBound* backCBound() {return _bounds[0].get();}
        
        virtual CMonomer* frontCMonomer() {return _monomers[_length - 1].get();}
        
        virtual CBound* frontCBound() {return _bounds[_length - 1].get();}
        
        
        ///Add a filament reaction to this subfilament
        virtual void addReaction(ReactionBase* r) {_reactions.push_back(r);}
        
        ///Get list of reactions associated with this subfilament
        virtual std::vector<ReactionBase*>& getReactions() {return _reactions;}
        
        ///Update all reactions associated with this subfilament
        virtual void updateReactions()
        {
            ///loop through all reactions, passivate/activate
            auto reactions = getReactions();
            
            for(auto &r : reactions) {
                
                if(r->getProductOfReactants() == 0)
                    r->passivateReaction();
                else
                    r->activateReaction();
            }
        }

        ///Get the current length
        virtual short length() {return _length;}
        
        ///Increase length
        virtual void increaseLength() {if(_length != _max_length) _length++;}
        
        ///Decrease length
        virtual void decreaseLength() {if(_length != 0) _length--;}
        
        ///see if the subfilament is at maxlength
        virtual bool atMaxLength() {return _length == _max_length;}
        
        ///set the length of this subfilament
        virtual void setLength(int length) {_length = length;}
        
        ///Print CSubFilament
        virtual void printCSubFilament()
        {
            std::cout << "Compartment:" << _compartment << std::endl;
            
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
        Filament* _mFilament;
        int _num_compartments = 1; ///number of compartments this filament spans
        int _length = 0; ///Length of filament
        
    public:
        ///Default constructor, does nothing
        CFilament() {};
        
        ///Default destructor, removes all Csubfilaments implicitly
        ~CFilament() {}
        
        ///Set mFilament ptr
        virtual void setMFilament(Filament* mFilament) { _mFilament = mFilament; }
        
        ///Add a CSubFilament
        virtual void addCSubFilament(CSubFilament* s) { addChild(std::unique_ptr<Component>(s)); }
        
        ///Remove a CSubFilament
        virtual void removeCSubFilament(CSubFilament* s) { removeChild(s); }
        
        ///Get front subfilament unconditionally (could be empty)
        ///@note -  no check on the number of children
        virtual CSubFilament* getFrontCSubFilament()
        {
            return static_cast<CSubFilament*>(children(numberOfChildren() - 1));
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
        ///@note should only be used when initializing filament
        virtual void setLength(int length) {
            _length = length;
        }
        
        ///get length of filament
        virtual int length() {return _length;}
        
        ///increase the number of compartments
        virtual void increaseNumCompartments() {_num_compartments++;}

        ///Get the number of compartments this filament spans
        virtual int numCompartments() {return _num_compartments;}
        
        ///Update all reactions associated with subfilament children
        virtual void updateReactions()
        {
            for(auto &c : children())
                static_cast<CSubFilament*>(c.get())->updateReactions();
        }

        ///Print entire filament
        virtual void printCFilament() {
            
            std::cout << "Filament ptr:" << this << std::endl;
            std::cout << "Length = " << _length <<std::endl;
            int index = 0;
            for (auto &c : children()) {
                std::cout << "CSubFilament " << index++ << ":" <<std::endl;
                static_cast<CSubFilament*>(c.get())->printCSubFilament();
                std::cout<< std::endl<<std::endl;
            }
        }
        
    };
    

    
}; //namespace chem



#endif /* defined(__CytoSim__CFilament__) */
