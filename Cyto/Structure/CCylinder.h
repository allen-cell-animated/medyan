//
//  CCylinder.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CCylinder__
#define __CytoSim__CCylinder__

#include <iostream>
#include "Compartment.h"
#include "CFilamentElement.h"

class ChemSimReactionKey;
class Cylinder;

/// CCylinder class holds all Monomers and Bounds
/*! 
 *  The CCylinder Class is an template class that has lists of the monomers and bounds that it contains.
 *  it has functionality to print the current composition.
 *  Accessing a particular species in the CCylinder is possible as well.
 */
class CCylinder {
    
protected:
    std::vector<std::unique_ptr<CMonomer>> _monomers; ///< list of monomers in this ccylinder
    std::vector<std::unique_ptr<CBound>> _bounds; ///< list of bound species in this ccylinder
    std::vector<ReactionBase*> _reactions;///< list of reactions associated with ccylinder
    Compartment* _compartment; ///< compartment this ccylinder is in
    
    Cylinder* _pCylinder;

public:
    ///Default constructor, sets compartment
    CCylinder(Compartment* c) : _compartment(c) {}
    
    /// Copy constructor
    /// @note This constructor will create a new ccylinder with different species and reactions
    /// within the compartment that is chosen as a parameter to the constructor. The copied and
    /// original ccylinders will not share reactions or species, but will be copied into a new compartment.
    CCylinder(const CCylinder& rhs, Compartment* c) : _compartment(c)
    {
        ///copy all monomers, bounds
        for(auto &m : rhs._monomers)
            _monomers.push_back(std::unique_ptr<CMonomer>(m->clone(c)));
        for(auto &b : rhs._bounds)
            _bounds.push_back(std::unique_ptr<CBound>(b->clone(c)));
        
        ///copy all reactions
        for(auto &r: rhs._reactions) {
            ReactionBase *R = _compartment->addInternalReactionUnique(
                                std::unique_ptr<ReactionBase>(r->clone(c->speciesContainer())));
            _reactions.push_back(R);
        }
    }
    
    /// Assignment is not allowed
    CCylinder& operator=(CCylinder &rhs) = delete;
    
    ///Default destructor, explicitly removes monomers and bounds (including their species, rxns)
    ~CCylinder()
    {
        ///Remove all reactions
        for(auto &r: _reactions)
            _compartment->removeInternalReaction(r);
        
        ///Remove all species
        for(auto &m: _monomers)
            for(auto &s : m->species())
                _compartment->removeSpecies(s);

        for(auto &b : _bounds)
            for(auto &s : b->species())
                _compartment->removeSpecies(s);
        
        _monomers.clear();
        _bounds.clear();
    }
    
    ///Clone, calls copy constructor
    virtual CCylinder* clone(Compartment* c) {
        return new CCylinder(*this, c);
    }
    
    ///get filament compartment
    Compartment* compartment() {return _compartment;}
    
    ///set parent cylinder
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    ///get parent cylinder
    Cylinder* getCylinder() {return _pCylinder;}
    
    ///Add a monomer to this CCylinder
    virtual void addCMonomer(CMonomer* monomer) { _monomers.emplace_back(std::unique_ptr<CMonomer>(monomer));}
    ///Add a bound to this CCylinder
    virtual void addCBound(CBound* bound) {_bounds.emplace_back(std::unique_ptr<CBound>(bound));}
    
    ///Get monomer at an index
    ///@note no check on index
    virtual CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    ///Get bound at an index
    ///@note no check on index
    virtual CBound* getCBound(int index) {return _bounds[index].get();}
    
    ///Add a filament reaction to this CCylinder
    virtual void addReaction(ReactionBase* r) {_reactions.push_back(r);}
    ///Get list of reactions associated with this CCylinder
    virtual std::vector<ReactionBase*>& getReactions() {return _reactions;}
    
    ///Add all reactions associated with this CCylinder
    virtual void addChemSimReactions()
    {
        for (auto &r: _reactions)
            ChemSim::addReaction(ChemSimReactionKey(), r);
    }
    
    ///Update all reactions associated with this CCylinder
    virtual void updateReactions()
    {
        ///loop through all reactions, passivate/activate
        for(auto &r : _reactions) {
            
            if(r->getProductOfReactants() == 0)
                r->passivateReaction();
            else
                r->activateReaction();
        }
    }
    
    ///Print CCylinder
    virtual void printCCylinder()
    {
        std::cout << "Compartment:" << _compartment << std::endl;
        
        std::cout << "Composition of CCylinder: " << std::endl;
        for (auto &m : _monomers){
            m->print();
            std::cout << ":";
        }
        std::cout << std::endl << "Bounds of CCylinder: " <<std::endl;
        for (auto &b : _bounds) {
            b->print();
            std::cout << ":";
        }
    }
};



#endif /* defined(__CytoSim__CCylinder__) */
