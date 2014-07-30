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
#include "CompartmentContainer.h"

namespace chem {
    class CMonomer;
    class CBound;
    
    /// CCylinder class holds all Monomers and Bounds
    /*! 
     *  The CCylinder Class is an template class that has lists of the monomers and bounds that it contains.
     *  it has functionality to print the current composition.
     *  Accessing a particular species in the CCylinder is possible as well.
     */
    class CCylinder : public Component {
        
    protected:
        std::vector<std::unique_ptr<CMonomer>> _monomers; ///< list of monomers in this sub filament
        std::vector<std::unique_ptr<CBound>> _bounds; ///< list of bound species in this sub filament
        std::vector<ReactionBase*> _reactions;///< list of reactions associated with this subfilament
        Compartment* _compartment; ///< compartment this CCylinder is in
        
    public:
        ///Default constructor, sets compartment
        CCylinder(Compartment* c) : _compartment(c) {}
        
        ///Default destructor, explicitly removes monomers and bounds
        ~CCylinder()
        {
            _monomers.clear();
            _bounds.clear();
        }
        
        ///get filament compartment
        Compartment* compartment() {return _compartment;}
        
        ///Add a monomer to this CCylinder
        virtual void addCMonomer(CMonomer* monomer) {
            _monomers.emplace_back(std::unique_ptr<CMonomer>(monomer));
        }
        
        ///Add a bound to this CCylinder
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
        
        virtual CMonomer* frontCMonomer() {return _monomers.back().get();}
        
        virtual CBound* frontCBound() {return _bounds.back().get();}
        
        ///Add a filament reaction to this subfilament
        virtual void addReaction(ReactionBase* r) {_reactions.push_back(r);}
        
        ///Get list of reactions associated with this subfilament
        virtual std::vector<ReactionBase*>& getReactions() {return _reactions;}
        
        ///Update all reactions associated with this subfilament
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
    
    
}; //namespace chem



#endif /* defined(__CytoSim__CCylinder__) */
