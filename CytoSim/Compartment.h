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
#include <unordered_map>
#include "Reaction.h"
#include "Species.h"
#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "Composite.h"

namespace chem {
    
    class Compartment : public Composite {
    protected:
        SpeciesPtrContainerVector _species;
        ReactionPtrContainerVector _internal_reactions;
        ReactionPtrContainerVector _diffusion_reactions;
        std::vector<Compartment*> _neighbours;
        std::unordered_map<int,float> _diffusion_rates;
    public:
        Compartment() : _species(), _internal_reactions(), _diffusion_reactions(), _neighbours(), _diffusion_rates() {}
        
        Compartment(const Compartment &C) : _species(), _internal_reactions(), _diffusion_reactions(), _neighbours(), _diffusion_rates()
        {
            C.cloneSpecies(this);
            C.cloneReactions(this);
            _diffusion_rates = C._diffusion_rates;
        }
        
//        friend void swap(Compartment& first, Compartment& second) // nothrow
//        {
//            // enable ADL (not necessary in our case, but good practice)
//            using std::swap;
//            swap(first._species, second._species);
//            for(auto &s : first._species.species()) s->setParent(&first);
//            for(auto &s : second._species.species()) s->setParent(&second);
//            swap(first._internal_reactions, second._internal_reactions);
//            for(auto &r : first._internal_reactions.reactions()) r->setParent(&first);
//            for(auto &r : second._internal_reactions.reactions()) r->setParent(&second);
//            swap(first._diffusion_reactions, second._diffusion_reactions);
//            for(auto &r : first._diffusion_reactions.reactions()) r->setParent(&first);
//            for(auto &r : second._diffusion_reactions.reactions()) r->setParent(&second);
//            swap(first._neighbours, second._neighbours);
//            swap(first._diffusion_rates, second._diffusion_rates);
//        }
        
        Compartment& operator=(const Compartment &other)
        {
//            swap(*this, other);
            _species.clear();
            _internal_reactions.clear();
            _diffusion_reactions.clear();
            other.cloneSpecies(this);
            other.cloneReactions(this);
            _diffusion_rates = other._diffusion_rates;
            return *this;
            // Note that _neighbours is not copied
        }

        
        virtual ~Compartment() noexcept {
            _internal_reactions.reactions().clear();
            _diffusion_reactions.reactions().clear();
            _species.species().clear();
        }
        
        virtual bool isSpeciesContainer() const {return true;}
        virtual bool isReactionsContainer() const {return true;}
        virtual std::string getFullName() const {return "Compartment";};
        size_t numberOfSpecies() const {return _species.species().size();}
        size_t numberOfReactions() const {return _internal_reactions.reactions().size()+_diffusion_reactions.reactions().size();}
        
        virtual Species* findSpeciesByName(const std::string &name) {return _species.findSpeciesByName(name);};
        virtual Species* findSpeciesByIndex (size_t index) {return _species.findSpeciesByIndex(index);};
        virtual Species* findSpeciesByMolecule (int molecule) {return _species.findSpeciesByMolecule(molecule);};
        
        virtual Species* addSpeciesUnique (std::unique_ptr<Species> &&species, float diff_rate = -1.0)
        {
            Species *sp = _species.addSpeciesUnique(std::move(species));
            sp->setParent(this);
            _diffusion_rates[sp->getMolecule()]=diff_rate;
            return sp;
        }
        
        virtual Reaction* addInternalReactionUnique (std::unique_ptr<Reaction> &&reaction)
        {
            Reaction *r = _internal_reactions.addReactionUnique(std::move(reaction));
            r->setParent(this);
            return r;
        }
        
        virtual Reaction* addDiffusionReactionUnique (std::unique_ptr<Reaction> &&reaction)
        {
            Reaction *r = _diffusion_reactions.addReactionUnique(std::move(reaction));
            r->setParent(this);
            return r;
        }
        
        template<typename ...Args>
        Species* addSpecies(Args&& ...args)
        {
//            std::cout << "Compartment::addSpecies()..." << std::endl;
            Species *sp = _species.addSpecies<SpeciesDiffusing>(std::forward<Args>(args)...);
            sp->setParent(this);
            _diffusion_rates[sp->getMolecule()]=-1.0; // not clear yet how to add diff_rate as an argument to this method
            return sp;
        }
        
        template<typename ...Args>
        Reaction* addInternalReaction (Args&& ...args)
        {
            //            std::cout << "Compartment::addReaction()..." << std::endl;
            Reaction *r = _internal_reactions.addReaction(std::forward<Args>(args)...);
            r->setParent(this);
            return r;
        }
        
        template<typename ...Args>
        Reaction* addDiffusionReaction (Args&& ...args)
        {
            //            std::cout << "Compartment::addReaction()..." << std::endl;
            Reaction *r = _diffusion_reactions.addReaction(std::forward<Args>(args)...);
            r->setParent(this);
            return r;
        }
        
        void setDiffusionRate(Species *sp, float diff_rate)
        {
            int molecule = sp->getMolecule();
            _diffusion_rates[molecule]=diff_rate;
        }
        
        void setDiffusionRate(int molecule, float diff_rate)
        {
            _diffusion_rates[molecule]=diff_rate;
        }
        
        void setDiffusionRate(std::string species_name, float diff_rate)
        {
            int molecule = SpeciesNamesDB::Instance()->stringToInt(species_name);
            _diffusion_rates[molecule]=diff_rate;
        }

        void addNeighbour(Compartment *comp)
        {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit==_neighbours.end())
                _neighbours.push_back(comp);
            else
                throw std::runtime_error("Compartment::addNeighbour(): Compartment is already a neighbour");
        }
        
        void removeNeighbour(Compartment *comp)
        {
            auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
            if(nit!=_neighbours.end())
                _neighbours.erase(nit);
            else
                throw std::out_of_range("Compartment::removeNeighbour(): Compartment is not a neighbour");
        }
        
        void cloneSpecies (Compartment *target) const
        {
            assert(target->numberOfSpecies()==0);
            for(auto &s : _species.species()){
                target->addSpeciesUnique(std::unique_ptr<Species>(s->clone()));
            }
        }
        
        void cloneReactions (Compartment *target) const
        {
            assert(target->numberOfReactions()==0);
            for(auto &r : _internal_reactions.reactions()){
                target->addInternalReactionUnique(std::unique_ptr<Reaction>(r->clone(target->_species.species().begin(),target->_species.species().end())));
            }
        }

        void cloneSpeciesReactions(Compartment* C)
        {
            this->cloneSpecies(C);
            this->cloneReactions(C);
            C->_diffusion_rates = this->_diffusion_rates;
            
        }
        
        virtual Compartment* clone()
        {
            Compartment *C = new Compartment();
            cloneSpeciesReactions(C);
            return C;
        }
        
        void generateDiffusionReactions()
        {
            for(auto &sp_this : _species.species()) {
                int molecule = sp_this->getMolecule();
                int diff_rate = _diffusion_rates[molecule];
                if(diff_rate<0) // Based on a convention that diffusing reactions require positive rates
                    continue;
                for(auto &C : _neighbours){
                    Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
                    Reaction *R = new Reaction({sp_this.get(),sp_neighbour},1,1,diff_rate);
                    this->addDiffusionReactionUnique(std::unique_ptr<Reaction>(R));
                }
            }
        }
        
        virtual Species* findSimilarSpecies (const Species &s) {return _species.findSimilarSpecies(s);}

        virtual Reaction* findSimilarInternalReaction (const Reaction &r)
        {
            return _internal_reactions.findSimilarReaction(r);
        }
        
        virtual Reaction* findSimilarDiffusionReaction (const Reaction &r)
        {
            return _diffusion_reactions.findSimilarReaction(r);
        }
            
        size_t numberOfNeighbours() const {return _neighbours.size();}
        
        std::vector<Compartment*>& neighbours() {return _neighbours;}
        const std::vector<Compartment*>& neighbours() const {return _neighbours;}
        
        virtual void printSpecies() {_species.printSpecies();}
        virtual void printReactions()
        {
            _internal_reactions.printReactions();
            _diffusion_reactions.printReactions();
        }
        
        virtual bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}

    };



}// end of chem
//

#endif
