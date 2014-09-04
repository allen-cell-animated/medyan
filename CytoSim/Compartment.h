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

class Bead;
class BoundaryElement;
class ChemSimReactionKey;


/// Compartment class is a compartment for holding species and reactions

/*! The Compartment class is a container for species, internal reactions, and diffusion
 *  reactions that can occur. A compartment object keeps track of the above while also
 *  holding pointers to its neighbors in order to generate diffusion reactions and other
 *  cross-compartment reactions.
 *
 *  Compartment initialization looks like the following:
 *  @code
 *  Compartment *C1 = new Compartment();
 *  Species *A = C1->addSpecies("A",99U);
 *  C1->setDiffusionRate(A,2000);
 *  @endcode
 */

class Compartment : public Composite {
protected:
    ///Reaction-diffusion components
    SpeciesPtrContainerVector _species;  ///< Container with all species in this compartment
    ReactionPtrContainerVector _internal_reactions; ///< Container with all internal reactions in compartment
    ReactionPtrContainerVector _diffusion_reactions; ///< Container with all diffusion reactions in compartment
    std::vector<Compartment*> _neighbours; ///< Neighbors of the compartment
    std::unordered_map<int,float> _diffusion_rates; ///< Diffusion rates of species in compartment
    
    std::vector<Bead*> _beads; ///<vector of beads that are in this compartment
    std::vector<BoundaryElement*> _boundaryElements; ///< vector of boundary element that are in this compartment

    bool _activated = false; ///< the compartment is activated for diffusion
    
public:
    ///Default constructor, only takes in number of dimensions
    Compartment() : _species(), _internal_reactions(), _diffusion_reactions(),
                    _neighbours(), _diffusion_rates() {}
    
    ///Constructor which clones another compartment
    Compartment(const Compartment &C) : _species(), _internal_reactions(), _diffusion_reactions(), _neighbours()
    {
        C.cloneSpecies(this);
        C.cloneReactions(this);
        _diffusion_rates = C._diffusion_rates;
        _activated = C._activated;
    }
    
    //Assignment operator
    Compartment& operator=(const Compartment &other);
    
    //Destructor
    virtual ~Compartment() noexcept
    {
        clearNeighbouringReactions();
        clearReactions();
        clearSpecies();
        removeFromNeighboursList();
    }
    
    /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain Species.
    virtual bool apply_impl(SpeciesVisitor &v) override;
    
    /// Applies ReactionVisitor v to every ReactionBase* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain ReactionBase.
    virtual bool apply_impl(ReactionVisitor &v) override;
    
    ///Activate a compartment
    virtual void activate() {_activated = true;}
    
    ///Deactivate a compartment
    virtual void deactivate() {_activated = false;}
    
    ///Check if compartment is activated
    virtual bool isActivated() {return _activated;}
    
    /// Removes all reactions from this compartment, diffusing and internal
    virtual void clearReactions()
    {
        _internal_reactions.reactions().clear();
        _diffusion_reactions.reactions().clear();
    }
    
    /// Clear all species from this compartment
    virtual void clearSpecies()
    {
        _species.species().clear();
    }

    /// Remove diffusion reactions between this compartment and its neighbors
    virtual void clearNeighbouringReactions()
    {
        for(auto &s : _species.species())
        {
            for(auto &n : _neighbours)
            {
                n->removeDiffusionReactions(s.get());
                n->removeInternalReactions(s.get());
            }
        }
    }
    /// Remove this compartment from the neighbor list
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
    
    /// This is a species container
    virtual bool isSpeciesContainer() const {return true;}
    /// This is a reaction container
    virtual bool isReactionsContainer() const {return true;}
    /// Returns compartment name
    virtual std::string getFullName() const {return "Compartment";};
    /// Returns the number of species in this compartment
    size_t numberOfSpecies() const {return _species.species().size();}
    /// Returns the number of internal reactions in this compartment
    size_t numberOfInternalReactions() const {return _internal_reactions.reactions().size();}
    /// Returns the total number of reactions in this compartment, diffusing and internal
    size_t numberOfReactions() const {return _internal_reactions.reactions().size()+_diffusion_reactions.reactions().size();}
    
    /// Species finder functions
    Species* findSpeciesByName(const std::string &name) {return _species.findSpeciesByName(name);};
    Species* findSpeciesDiffusingByName(const std::string &name) {return _species.findSpeciesDiffusingByName(name);};
    Species* findSpeciesByIndex (size_t index) {return _species.findSpeciesByIndex(index);};
    Species* findSpeciesByMolecule (int molecule) {return _species.findSpeciesByMolecule(molecule);};
    Species* findSimilarSpecies (const Species &s) {return _species.findSimilarSpecies(s);}
    
    ///Remove species from this compartment
    size_t removeSpecies(Species* species) {return _species.removeSpecies(species);}
    
    /// Finds a similar internal reaction, see ReactionBase function
    ReactionBase* findSimilarInternalReaction (const ReactionBase &r)
    {
        return _internal_reactions.findSimilarReaction(r);
    }
    
    /// Finds a similar diffusion reaction, see ReactionBase function
    ReactionBase* findSimilarDiffusionReaction (const ReactionBase &r)
    {
        return _diffusion_reactions.findSimilarReaction(r);
    }
    
    /// Remove all diffusion reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeDiffusionReactions (Species* s) {_diffusion_reactions.removeReactions(s);}
    
    /// Remove all internal reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeInternalReactions (Species* s) {_internal_reactions.removeReactions(s);}
    
    /// Remove a diffusion reaction
    virtual void removeDiffusionReaction(ReactionBase *r) {_diffusion_reactions.removeReaction(r);}
    
    /// Remove an internal reaction
    virtual void removeInternalReaction(ReactionBase *r) {_internal_reactions.removeReaction(r);}
    
    /// Add a unique species pointer to this compartment
    Species* addSpeciesUnique (std::unique_ptr<Species> &&species, float diff_rate = -1.0)
    {
        Species *sp = _species.addSpeciesUnique(std::move(species));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=diff_rate;
        return sp;
    }
    
    /// Add a unique internal reaction pointer to this compartment
    ReactionBase* addInternalReactionUnique (std::unique_ptr<ReactionBase> &&reaction)
    {
        ReactionBase *r = _internal_reactions.addReactionUnique(std::move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Add a unique diffusing reaction pointer to this compartment
    ReactionBase* addDiffusionReactionUnique (std::unique_ptr<ReactionBase> &&reaction)
    {
        ReactionBase *r = _diffusion_reactions.addReactionUnique(std::move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Add a diffusing species to this compartment
    /// @param args - any number of SpeciesDiffusing objects
    template<typename ...Args>
    Species* addSpecies(Args&& ...args)
    {
        Species *sp = _species.addSpecies<SpeciesDiffusing>(std::forward<Args>(args)...);
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0; // not clear yet how to add diff_rate as an argument to this method
        return sp;
    }
    
    /// Add a filament species to this compartment
    /// @param args - any number of SpeciesDiffusing objects
    template<typename ...Args>
    SpeciesFilament* addSpeciesFilament(Args&& ...args)
    {
        SpeciesFilament *sp =
            static_cast<SpeciesFilament*>(_species.addSpecies<SpeciesFilament>(std::forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0; // not clear yet how to add diff_rate as an argument to this method
        return sp;
    }
    
    /// Add a bound species to this compartment
    /// @param args - any number of SpeciesDiffusing objects
    template<typename ...Args>
    SpeciesBound* addSpeciesBound(Args&& ...args)
    {
        SpeciesBound *sp =
            static_cast<SpeciesBound*>(_species.addSpecies<SpeciesBound>(std::forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0; // not clear yet how to add diff_rate as an argument to this method
        return sp;
    }
    

    /// Add an internal reaction to this compartment
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addInternalReaction (Args&& ...args)
    {
        //            std::cout << "Compartment::addReaction()..." << std::endl;
        ReactionBase *r = _internal_reactions.addReaction<M,N>(std::forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    /// Add an internal reaction to this compartment
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addInternal(std::initializer_list<Species*> species, float rate)
    {
        ReactionBase *r = _internal_reactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }
    
    /// Add a diffusion reaction to this compartment
    template<typename ...Args>
    ReactionBase* addDiffusionReaction (Args&& ...args)
    {
        //            std::cout << "Compartment::addReaction()..." << std::endl;
        ReactionBase *r = _diffusion_reactions.addReaction<1,1>(std::forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    
    ///Add a bead to this compartment
    void addBead(Bead* b) {_beads.push_back(b);}
    
    ///Remove a bead from this compartment
    ///@note does nothing if bead is not in compartment already
    void removeBead(Bead* b) {
        auto it = std::find(_beads.begin(), _beads.end(), b);
        if(it != _beads.end()) _beads.erase(it);
    }
    
    ///get the beads in this compartment
    std::vector<Bead*>& getBeads() {return _beads;}
    
    ///Add a boundary element to this compartment
    void addBoundaryElement(BoundaryElement* be) {_boundaryElements.push_back(be);}
    
    ///Remove a boundary element from this compartment
    ///@note does nothing if boundary element is not in compartment
    void removeBoundaryElement(BoundaryElement* be) {
        auto it = std::find(_boundaryElements.begin(), _boundaryElements.end(), be);
        if(it != _boundaryElements.end()) _boundaryElements.erase(it);
    }
    
    ///get the boundary elements in this compartment
    std::vector<BoundaryElement*>& getBoundaryElements() {return _boundaryElements;}
    
    /// Get the diffusion rate of a species
    /// @param - species_name, a string
    float getDiffusionRate(std::string species_name)
    {
        int molecule = SpeciesNamesDB::Instance()->stringToInt(species_name);
        return _diffusion_rates[molecule];
    }
    
    /// Set the diffusion rate of a species in the compartment
    void setDiffusionRate(Species *sp, float diff_rate)
    {
        int molecule = sp->getMolecule();
        _diffusion_rates[molecule]=diff_rate;
    }
    
    /// Set the diffusion rate of a species in the compartment
    /// @param - molecule, an integer representing the species in the speciesDB
    void setDiffusionRate(int molecule, float diff_rate)
    {
        _diffusion_rates[molecule]=diff_rate;
    }
    
    /// Set the diffusion rate of a species in the compartment
    /// @param - species_name, a string
    void setDiffusionRate(std::string species_name, float diff_rate)
    {
        int molecule = SpeciesNamesDB::Instance()->stringToInt(species_name);
        _diffusion_rates[molecule]=diff_rate;
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void addNeighbour(Compartment *comp)
    {
        auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit==_neighbours.end())
            _neighbours.push_back(comp);
        else
            throw std::runtime_error("Compartment::addNeighbour(): Compartment is already a neighbour");
    }
    
    /// Remove a neighboring compartment
    void removeNeighbour(Compartment *comp)
    {
        auto nit = std::find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit!=_neighbours.end())
            _neighbours.erase(nit);
        else
            throw std::out_of_range("Compartment::removeNeighbour(): Compartment is not a neighbour");
    }
    
    /// Clone the species values of another compartment into this one
    void cloneSpecies (Compartment *target) const
    {
        assert(target->numberOfSpecies()==0);
        for(auto &s : _species.species()){
            target->addSpeciesUnique(std::unique_ptr<Species>(s->clone()));
        }
    }
    
    /// Clone the reaction values of another compartment into this one
    void cloneReactions (Compartment *target) const
    {
        assert(target->numberOfReactions()==0);
        for(auto &r : _internal_reactions.reactions()){
            target->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r->clone(target->_species)));
        }
    }

    /// Clone both the species and compartments into this compartment
    void cloneSpeciesReactions(Compartment* C)
    {
        this->cloneSpecies(C);
        this->cloneReactions(C);
        C->_diffusion_rates = this->_diffusion_rates;
        
    }
    
    /// Clone a compartment
    /// @note - this does not clone the neighbors, just reactions and species
    virtual Compartment* clone()
    {
        Compartment *C = new Compartment(*this);
        return C;
    }
    
    /// Generate diffusion reactions between this compartment and another
    ///@return - a vector of reactionbases that was just added 
    std::vector<ReactionBase*> generateDiffusionReactions(Compartment* C);

    /// Generate all diffusion reactions for this compartment and its neighbors
    void generateAllDiffusionReactions();
    
    /// Gives the number of neighbors to this compartment
    size_t numberOfNeighbours() const {return _neighbours.size();}
    
    
    ///Get the species container vector
    SpeciesPtrContainerVector& speciesContainer() {return _species;}
    const SpeciesPtrContainerVector& speciesContainer() const {return _species;}
    
    ///Get the internal reaction container vector
    ReactionPtrContainerVector& internalReactionContainer() {return _internal_reactions;}
    const ReactionPtrContainerVector& internalReactionContainer() const {return _internal_reactions;}
    
    ///Get the diffusion reaction container vector
    ReactionPtrContainerVector& diffusionReactionContainer() {return _diffusion_reactions;}
    const ReactionPtrContainerVector& diffusionReactionContainer() const {return _diffusion_reactions;}
    
    /// Get the vector list of neighbors to this compartment
    std::vector<Compartment*>& neighbours() {return _neighbours;}
    const std::vector<Compartment*>& neighbours() const {return _neighbours;}
    
    /// Print the species in this compartment
    void printSpecies() {_species.printSpecies();}
    /// Print the reactions in this compartment
    void printReactions()
    {
        _internal_reactions.printReactions();
        _diffusion_reactions.printReactions();
    }
    
    /// Check if all species are unique pointers
    bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}
    
    /// Adds the reactions of this compartment to the ChemSim object
    /// @param - chem, a ChemSim object that runs the reaction-diffusion algorithm
    virtual void addChemSimReactions()
    {
        for(auto &r1 : _internal_reactions.reactions())
            ChemSim::addReaction(ChemSimReactionKey(), r1.get());
        for(auto &r2 : _diffusion_reactions.reactions())
            ChemSim::addReaction(ChemSimReactionKey(), r2.get());
    }
    
    /// Print properties of this compartment
    virtual void printSelf()
    {
        std::cout << this->getFullName() << "\n"
        << "Number of neighbors: " << numberOfNeighbours() << std::endl;
        printSpecies();
        std::cout << "Reactions:" << std::endl;
        printReactions();
    }

};
#endif
