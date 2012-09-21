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
#include <array>
#include <unordered_map>
#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "Composite.h"
#include "ChemSim.h"

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
        
        Compartment& operator=(const Compartment &other);
        
        virtual ~Compartment() noexcept
        {
//            std::cout << "~Compartment(): ptr=" << this << std::endl;
            clearReactions();
            clearNeighbouringReactions();
            clearSpecies();
            removeFromNeighboursList();
        }

        /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
        /// This method needs to be overriden by descendent classes that contain Species.
        virtual bool apply_impl(SpeciesVisitor &v) override;
        
        /// Applies ReactionVisitor v to every ReactionBase* object directly owned by this node.
        /// This method needs to be overriden by descendent classes that contain ReactionBase.
        virtual bool apply_impl(ReactionVisitor &v) override;
        
        virtual void clearReactions()
        {
            _internal_reactions.reactions().clear();
            _diffusion_reactions.reactions().clear();
        }
        
        virtual void clearSpecies()
        {
            _species.species().clear();
        }
        
        virtual void clearNeighbouringReactions()
        {
            for(auto &s : _species.species())
            {
                for(auto &n : _neighbours)
                {
                    n->removeDiffusionReactions(s.get());
                }
            }
        }
            
        virtual void removeFromNeighboursList()
        {
            for(auto &n : _neighbours)
            {
                n->removeNeighbour(this);
            }
        }

        
        /// Returns true if two Compartment objects are equal.
        /// Two Compartment objects are equal if each contains analogous Species and Reaction objects, in the same order
        friend bool operator==(const Compartment& a, const Compartment& b);
        
        /// Return true if two Compartment are not equal.
        /// @see operator ==(const Compartment& a, const Compartment& b) above
        friend bool operator !=(const Compartment& a, const Compartment& b)
        {
            return !(a==b);
        }
        
        virtual bool isSpeciesContainer() const {return true;}
        virtual bool isReactionsContainer() const {return true;}
        virtual std::string getFullName() const {return "Compartment";};
        size_t numberOfSpecies() const {return _species.species().size();}
        size_t numberOfInternalReactions() const {return _internal_reactions.reactions().size();}
        size_t numberOfReactions() const {return _internal_reactions.reactions().size()+_diffusion_reactions.reactions().size();}
        
        Species* findSpeciesByName(const std::string &name) {return _species.findSpeciesByName(name);};
        Species* findSpeciesByIndex (size_t index) {return _species.findSpeciesByIndex(index);};
        Species* findSpeciesByMolecule (int molecule) {return _species.findSpeciesByMolecule(molecule);};
        Species* findSimilarSpecies (const Species &s) {return _species.findSimilarSpecies(s);}
        ReactionBase* findSimilarInternalReaction (const ReactionBase &r)
        {
            return _internal_reactions.findSimilarReaction(r);
        }
        ReactionBase* findSimilarDiffusionReaction (const ReactionBase &r)
        {
            return _diffusion_reactions.findSimilarReaction(r);
        }
        
        virtual void removeDiffusionReactions (Species* s) {_diffusion_reactions.removeReactions(s);}
        
        Species* addSpeciesUnique (std::unique_ptr<Species> &&species, float diff_rate = -1.0)
        {
            Species *sp = _species.addSpeciesUnique(std::move(species));
            sp->setParent(this);
            _diffusion_rates[sp->getMolecule()]=diff_rate;
            return sp;
        }
        
        ReactionBase* addInternalReactionUnique (std::unique_ptr<ReactionBase> &&reaction)
        {
            ReactionBase *r = _internal_reactions.addReactionUnique(std::move(reaction));
            r->setParent(this);
            return r;
        }
        
        ReactionBase* addDiffusionReactionUnique (std::unique_ptr<ReactionBase> &&reaction)
        {
            ReactionBase *r = _diffusion_reactions.addReactionUnique(std::move(reaction));
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
        
        template<unsigned short M, unsigned short N, typename ...Args>
        ReactionBase* addInternalReaction (Args&& ...args)
        {
            //            std::cout << "Compartment::addReaction()..." << std::endl;
            ReactionBase *r = _internal_reactions.addReaction<M,N>(std::forward<Args>(args)...);
            r->setParent(this);
            return r;
        }
        
        template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
        ReactionBase* addInternal(std::initializer_list<Species*> species, float rate)
        {
            ReactionBase *r = _internal_reactions.add<RXN,M,N>(species,rate);
            r->setParent(this);
            return r;
        }
        
        template<typename ...Args>
        ReactionBase* addDiffusionReaction (Args&& ...args)
        {
            //            std::cout << "Compartment::addReaction()..." << std::endl;
            ReactionBase *r = _diffusion_reactions.addReaction<1,1>(std::forward<Args>(args)...);
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
                target->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r->clone(target->_species)));
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
        
        void generateDiffusionReactions();
            
        size_t numberOfNeighbours() const {return _neighbours.size();}
        
        std::vector<Compartment*>& neighbours() {return _neighbours;}
        const std::vector<Compartment*>& neighbours() const {return _neighbours;}
        
        void printSpecies() {_species.printSpecies();}
        void printReactions()
        {
            _internal_reactions.printReactions();
            _diffusion_reactions.printReactions();
        }
        
        bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}
        
        virtual void addChemSimReactions(ChemSim &chem)
        {
            for(auto &r1 : _internal_reactions.reactions())
                chem.addReaction(r1.get());
            for(auto &r2 : _diffusion_reactions.reactions())
                chem.addReaction(r2.get());
        }
        
        virtual void printSelf()
        {
            std::cout << this->getFullName() << "\n"
            << "Number of neighbors: " << numberOfNeighbours() << std::endl;
            printSpecies();
            std::cout << "Reactions:" << std::endl;
            printReactions();
        }

    };

    template <size_t NDIM>
    class CompartmentCubic : public Compartment {
    private:
        std::array<float, NDIM> _coords;
        std::array<float, NDIM> _sides;
    public:
        CompartmentCubic() : Compartment() {};

        CompartmentCubic(const CompartmentCubic &other) : Compartment(other), _sides(other._sides) {}
        
        CompartmentCubic& operator=(const CompartmentCubic &other) {
            this->Compartment::operator=(other);
            _sides=other._sides;
        }
        
        virtual ~CompartmentCubic() noexcept {
        }

        virtual std::string getFullName() const {return "CompartmentCubic<"+std::to_string(NDIM)+">";};
        
        void setSideLength(size_t i, float value) {_sides[i]=value;}
        
        void setCoord(size_t i, float value) {_coords[i]=value;}

        template <class input_iter>
        void setSides(input_iter input_it)
        {
            std::copy_n(input_it,NDIM,_sides.begin());
        }
        
        template <class input_iter>
        void setCoords(input_iter input_it)
        {
            std::copy_n(input_it,NDIM,_coords.begin());
        }

        virtual void printSelf()
        {
            Compartment::printSelf();
            std::cout << "Coords: ";
            for(auto &x : _coords) std::cout << x << " ";
            std::cout << "\nSide: ";
            for(auto &y : _sides) std::cout << y << " ";
        }
        
        std::array<float, NDIM>& coords() {return _coords;}
        const std::array<float, NDIM>& coords() const {return _coords;}
        std::array<float, NDIM>& sides() {return _sides;}
        const std::array<float, NDIM>& sides() const {return _sides;}        
    };

}// end of chem
//

#endif
