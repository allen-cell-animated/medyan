//
//  Compartment.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Compartment_h
#define CytoSim_Experimenting_Compartment_h

#include <vector>
#include "Reaction.h"
#include "Species.h"
#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "Composite.h"

namespace chem {
    
    class Compartment : public Composite, public SpeciesPtrContainerVector, public ReactionPtrContainerVector {
    protected:
        std::vector<Compartment*> _neighbours;
    public:
        Compartment() :  Composite(), SpeciesPtrContainerVector() {}
        virtual ~Compartment() noexcept {
            _reactions.clear();
            _species.clear();
        }
        
        virtual bool isSpeciesContainer() const {return true;}
        virtual bool isReactionsContainer() const {return true;}
        virtual std::string getFullName() const {return "Compartment";};
        size_t numberOfSpecies() const {return _species.size();}
        size_t numberOfReactions() const {return _reactions.size();}
        
        virtual Species* addSpeciesUnique (std::unique_ptr<Species> &&species) {
            Species *sp = SpeciesPtrContainerVector::addSpeciesUnique(std::move(species));
            sp->setParent(this);
            return sp;
        }
        
        virtual Reaction* addReactionUnique (std::unique_ptr<Reaction> &&reaction) {
            Reaction *r = ReactionPtrContainerVector::addReactionUnique(std::move(reaction));
            r->setParent(this);
            return r;
        }
        
        template<typename ...Args>
        Species* addSpecies(Args&& ...args){
//            std::cout << "Compartment::addSpecies()..." << std::endl;
            Species *sp = SpeciesPtrContainerVector::addSpecies<SpeciesDiffusing>(std::forward<Args>(args)...);
            sp->setParent(this);
            return sp;
        }
        
        template<typename ...Args>
        Reaction* addReaction (Args&& ...args){
            //            std::cout << "Compartment::addReaction()..." << std::endl;
            Reaction *r = ReactionPtrContainerVector::addReaction(std::forward<Args>(args)...);
            r->setParent(this);
            return r;
        }
        
        std::vector<std::unique_ptr<Species>>& species() {return _species;}
        const std::vector<std::unique_ptr<Species>>& species() const {return _species;}

        std::vector<std::unique_ptr<Reaction>>& reactions() {return _reactions;}
        const std::vector<std::unique_ptr<Reaction>>& reactions() const {return _reactions;}

        void addNeighbour(Compartment *comp) {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit==_neighbours.end())
                _neighbours.push_back(comp);
            else
                throw std::runtime_error("Compartment::addNeighbour(): Compartment is already a neighbour");
        }
        
        void removeNeighbour(Compartment *comp) {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit!=_neighbours.end())
                _neighbours.erase(nit);
            else
                throw std::out_of_range("Compartment::removeNeighbour(): Compartment is not a neighbour");
        }
        
        void cloneSpecies (Compartment *target) {
            assert(target->species().size()==0);
            for(auto &s : _species){
                target->addSpeciesUnique(std::unique_ptr<Species>(s->clone()));
            }
        }
        
        void cloneReactions (Compartment *target) {
            assert(target->reactions().size()==0);
            for(auto &r : _reactions){
                target->addReactionUnique(std::unique_ptr<Reaction>(r->clone(target->species().begin(),target->species().end())));
            }
        }

        void cloneSpeciesReactions(Compartment* C) {
            this->cloneSpecies(C);
            this->cloneReactions(C);
        }
        
        virtual Compartment* clone() {
            Compartment *C = new Compartment();
            cloneSpeciesReactions(C);
            return C;
        }
        
        size_t numberOfNeighbours() const {return _neighbours.size();}
        
        std::vector<Compartment*>& neighbours() {return _neighbours;}
        
        const std::vector<Compartment*>& neighbours() const {return _neighbours;}
        
    };



}// end of chem
//

#endif
