
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Compartment_h
#define M3SYM_Compartment_h

#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>

#include "common.h"

#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "BindingManager.h"
#include "Composite.h"
#include "ChemSim.h"

//FORWARD DECLARATIONS
class Bead;
class Cylinder;
class BoundaryElement;

/// A container or holding Species and [Reactions](@ref Reactions).

/*! The Compartment class is a container for Species, [Reactions](@ref Reactions), and 
 *  diffusion [Reactions](@ref Reactions) that can occur. A Compartment object keeps 
 *  track of the above while also holding pointers to its neighbors in order to generate
 *  diffusion [Reactions](@ref Reactions) and other cross-compartment 
 *  [Reactions](@ref Reactions).
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
    
    ///CHEMICAL CONTAINERS
    SpeciesPtrContainerVector _species;  ///< Container with all species
                                         ///< in this compartment
    ReactionPtrContainerVector _internal_reactions; ///< Container with all internal
                                                    ///< reactions in compartment
    ReactionPtrContainerVector _diffusion_reactions; ///< Container with all diffusion
                                                     ///< reactions in compartment
    
    unordered_map<int,float> _diffusion_rates; ///< Diffusion rates of Species
                                               ///< in compartment
    
    /// All binding managers for this compartment
    vector<unique_ptr<FilamentBindingManager>> _bindingManagers;

    ///ELEMENT CONTAINERS
    unordered_set<BoundaryElement*> _boundaryElements; ///< Set of boundary element
                                                       ///< that are in this compartment
    
    unordered_set<Bead*> _beads; ///< Set of beads that are in this compartment
    
    unordered_set<Cylinder*> _cylinders; ///< Set of cylinders that are in this compartment
    
    vector<Compartment*> _neighbours; ///< Neighbors of the compartment

    
    ///OTHER COMPARTMENT PROPERTIES
    vector<double> _coords;  ///< Coordinates of this compartment
    bool _activated = false; ///< The compartment is activated for diffusion
    
public:
    /// Default constructor, only takes in number of dimensions
    Compartment() : _species(), _internal_reactions(), _diffusion_reactions(),
                    _neighbours(), _diffusion_rates() {}
    
    /// Constructor which clones another compartment
    Compartment(const Compartment &C) : _species(), _internal_reactions(),
                                        _diffusion_reactions(), _neighbours()
    {
        C.cloneSpecies(this);
        C.cloneReactions(this);
        _diffusion_rates = C._diffusion_rates;
        _activated = C._activated;
        // Should eventuall clone beads, cylinders, boundary elements.... not clear yet
    }
    
    /// Assignment operator
    Compartment& operator=(const Compartment &other);
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Compartment() noexcept
    {
        clearNeighbouringReactions();
        clearReactions();
        clearSpecies();
        removeFromNeighboursList();
        
        // Should eventually delete beads, cylinders, boundary elements....not yet clear
        
    }
    
    /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain Species.
    virtual bool apply_impl(SpeciesVisitor &v) override;
    
    /// Applies ReactionVisitor v to every ReactionBase* object directly owned by this
    /// node. This method needs to be overriden by descendent classes that contain
    /// ReactionBase.
    virtual bool apply_impl(ReactionVisitor &v) override;
    
    ///Activate a compartment
    virtual void activate() {_activated = true;}
    
    ///Deactivate a compartment
    virtual void deactivate() {_activated = false;}
    
    ///Check if compartment is activated
    virtual bool isActivated() {return _activated;}
    
    ///Setter and getter for coordinates
    virtual void setCoordinates(vector<double> coords) {_coords = coords;}
    virtual const vector<double>& coordinates() {return _coords;}
    
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
    /// Two Compartment objects are equal if each contains analogous Species and
    /// Reaction objects, in the same order
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
    virtual string getFullName() const {return "Compartment";};
    /// Returns the number of species in this compartment
    size_t numberOfSpecies() const {return _species.species().size();}
    /// Returns the number of internal reactions in this compartment
    size_t numberOfInternalReactions() const {
        return _internal_reactions.reactions().size();
    }
    /// Returns the total number of reactions in this compartment, diffusing and
    /// internal
    size_t numberOfReactions() const {
        return _internal_reactions.reactions().size() +
               _diffusion_reactions.reactions().size();
    }
    
    //@{
    /// Species finder functions
    Species* findSpeciesByName(const string &name) {
        return _species.findSpeciesByName(name);
    }
    Species* findSpeciesByIndex (size_t index) {
        return _species.findSpeciesByIndex(index);
    }
    Species* findSpeciesByMolecule (int molecule) {
        return _species.findSpeciesByMolecule(molecule);
    }
    Species* findSimilarSpecies (const Species &s) {
        return _species.findSimilarSpecies(s);
    }
    //@}
    
    ///Remove species from this compartment
    size_t removeSpecies(Species* species) {return _species.removeSpecies(species);}
    
    /// Finds a similar internal reaction, see ReactionBase function
    ReactionBase* findSimilarInternalReaction (const ReactionBase &r) {
        return _internal_reactions.findSimilarReaction(r);
    }
    
    /// Finds a similar diffusion reaction, see ReactionBase function
    ReactionBase* findSimilarDiffusionReaction (const ReactionBase &r) {
        return _diffusion_reactions.findSimilarReaction(r);
    }
    
    /// Remove all diffusion reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeDiffusionReactions (Species* s) {
        _diffusion_reactions.removeReactions(s);
    }
    
    /// Remove all internal reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeInternalReactions (Species* s) {
        _internal_reactions.removeReactions(s);
    }
    
    /// Remove a diffusion reaction
    virtual void removeDiffusionReaction(ReactionBase *r) {
        _diffusion_reactions.removeReaction(r);
    }
    
    /// Remove an internal reaction
    virtual void removeInternalReaction(ReactionBase *r) {
        _internal_reactions.removeReaction(r);
    }
    
    /// Add a unique species pointer to this compartment
    Species* addSpeciesUnique (unique_ptr<Species> &&species, float diff_rate = -1.0) {
        Species *sp = _species.addSpeciesUnique(move(species));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=diff_rate;
        return sp;
    }
    
    /// Add a unique internal reaction pointer to this compartment
    ReactionBase* addInternalReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _internal_reactions.addReactionUnique(move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Add a unique diffusing reaction pointer to this compartment
    ReactionBase* addDiffusionReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _diffusion_reactions.addReactionUnique(move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Add a diffusing species to this compartment
    /// @param args - any number of SpeciesDiffusing objects
    template<typename ...Args>
    SpeciesDiffusing* addSpeciesDiffusing(Args&& ...args) {
        SpeciesDiffusing *sp =
            static_cast<SpeciesDiffusing*>(
            _species.addSpecies<SpeciesDiffusing>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a filament species to this compartment
    /// @param args - any number of SpeciesFilament objects
    template<typename ...Args>
    SpeciesFilament* addSpeciesFilament(Args&& ...args) {
        SpeciesFilament *sp =
            static_cast<SpeciesFilament*>(
            _species.addSpecies<SpeciesFilament>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a plus end species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesPlusEnd* addSpeciesPlusEnd(Args&& ...args) {
        SpeciesPlusEnd *sp =
        static_cast<SpeciesPlusEnd*>(
        _species.addSpecies<SpeciesPlusEnd>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a minus end species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesMinusEnd* addSpeciesMinusEnd(Args&& ...args) {
        SpeciesMinusEnd *sp =
        static_cast<SpeciesMinusEnd*>(
        _species.addSpecies<SpeciesMinusEnd>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a bound species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesBound* addSpeciesBound(Args&& ...args) {
        SpeciesBound *sp =
        static_cast<SpeciesBound*>(
        _species.addSpecies<SpeciesBound>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a linker species to this compartment
    /// @param args - any number of SpeciesLinker objects
    template<typename ...Args>
    SpeciesLinker* addSpeciesLinker(Args&& ...args) {
        SpeciesLinker *sp =
        static_cast<SpeciesLinker*>(
            _species.addSpecies<SpeciesLinker>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a motor species to this compartment
    /// @param args - any number of SpeciesMotor objects
    template<typename ...Args>
    SpeciesMotor* addSpeciesMotor(Args&& ...args) {
        SpeciesMotor *sp =
        static_cast<SpeciesMotor*>(
            _species.addSpecies<SpeciesMotor>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a brancher species to this compartment
    /// @param args - any number of SpeciesBrancher objects
    template<typename ...Args>
    SpeciesBrancher* addSpeciesBrancher(Args&& ...args) {
        SpeciesBrancher *sp =
        static_cast<SpeciesBrancher*>(
            _species.addSpecies<SpeciesBrancher>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }
    
    /// Add a single binding species to this compartment
    /// @param args - any number of SpeciesSingleBinding objects
    template<typename ...Args>
    SpeciesSingleBinding* addSpeciesSingleBinding(Args&& ...args) {
        SpeciesSingleBinding *sb =
        static_cast<SpeciesSingleBinding*>(
            _species.addSpecies<SpeciesSingleBinding>(forward<Args>(args)...));
        sb->setParent(this);
        _diffusion_rates[sb->getMolecule()]=-1.0;
        return sb;
    }
    
    /// Add a pair binding species to this compartment
    /// @param args - any number of SpeciesPairBinding objects
    template<typename ...Args>
    SpeciesPairBinding* addSpeciesPairBinding(Args&& ...args) {
        SpeciesPairBinding *sb =
        static_cast<SpeciesPairBinding*>(
            _species.addSpecies<SpeciesPairBinding>(forward<Args>(args)...));
        sb->setParent(this);
        _diffusion_rates[sb->getMolecule()]=-1.0;
        return sb;
    }

    /// Add an internal reaction to this compartment
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addInternalReaction (Args&& ...args) {
        //            cout << "Compartment::addReaction()..." << endl;
        ReactionBase *r = _internal_reactions.addReaction<M,N>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    /// Add an internal reaction to this compartment
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addInternal(initializer_list<Species*> species, float rate) {
        ReactionBase *r = _internal_reactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }
    
    /// Add a diffusion reaction to this compartment
    template<typename ...Args>
    ReactionBase* addDiffusionReaction (Args&& ...args) {
        ReactionBase *r = _diffusion_reactions.addReaction<1,1>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    /// Add a binding manager to this compartment
    void addFilamentBindingManager(FilamentBindingManager* m) {
        _bindingManagers.emplace_back(m);
    }
    
    /// Get binding managers for this compartment
    vector<unique_ptr<FilamentBindingManager>>& getFilamentBindingManagers() {
        return _bindingManagers;
    }
    
    ///Add a bead to this compartment
    void addBead(Bead* b) {_beads.insert(b);}
    
    ///Remove a bead from this compartment
    ///@note does nothing if bead is not in compartment already
    void removeBead(Bead* b) {
        auto it = _beads.find(b);
        if(it != _beads.end()) _beads.erase(it);
    }
    
    ///get the beads in this compartment
   unordered_set<Bead*>& getBeads() {return _beads;}
    
    ///Add a boundary element to this compartment
    void addBoundaryElement(BoundaryElement* be) {_boundaryElements.insert(be);}
    
    ///Remove a boundary element from this compartment
    ///@note does nothing if boundary element is not in compartment
    void removeBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        if(it != _boundaryElements.end()) _boundaryElements.erase(it);
    }
    ///Check if boundary element is in this container
    bool hasBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        return (it != _boundaryElements.end());   
    }
    ///get the boundary elements in this compartment
   unordered_set<BoundaryElement*>& getBoundaryElements() {return _boundaryElements;}
    
    ///Add a cylinder to this compartment
    void addCylinder(Cylinder* c) {_cylinders.insert(c);}
    
    ///Remove a cylinder from this compartment
    ///@note does nothing if cylinder is not in compartment already
    void removeCylinder(Cylinder* c) {
        auto it = _cylinders.find(c);
        if(it != _cylinders.end()) _cylinders.erase(it);
    }
    ///get the cylinders in this compartment
   unordered_set<Cylinder*>& getCylinders() {return _cylinders;}
    
    /// Get the diffusion rate of a species
    /// @param - species_name, a string
    float getDiffusionRate(string species_name) {
        int molecule = SpeciesNamesDB::instance()->stringToInt(species_name);
        return _diffusion_rates[molecule];
    }
    
    /// Set the diffusion rate of a species in the compartment
    void setDiffusionRate(Species *sp, float diff_rate) {
        int molecule = sp->getMolecule();
        _diffusion_rates[molecule]=diff_rate;
    }
    
    /// Set the diffusion rate of a species in the compartment
    /// @param - molecule, an integer representing the species in the speciesDB
    void setDiffusionRate(int molecule, float diff_rate) {
        _diffusion_rates[molecule]=diff_rate;
    }
    
    /// Set the diffusion rate of a species in the compartment
    /// @param - species_name, a string
    void setDiffusionRate(string species_name, float diff_rate) {
        int molecule = SpeciesNamesDB::instance()->stringToInt(species_name);
        _diffusion_rates[molecule]=diff_rate;
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void addNeighbour(Compartment *comp) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit==_neighbours.end())
            _neighbours.push_back(comp);
        else
            throw runtime_error(
            "Compartment::addNeighbour(): Compartment is already a neighbour");
    }
    
    /// Remove a neighboring compartment
    void removeNeighbour(Compartment *comp) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit!=_neighbours.end())
            _neighbours.erase(nit);
    }
    
    /// Clone the species values of another compartment into this one
    void cloneSpecies (Compartment *target) const {
        assert(target->numberOfSpecies()==0);
        for(auto &s : _species.species()){
            Species* sClone = s->clone();
            target->addSpeciesUnique(unique_ptr<Species>(sClone));
        }
    }
    
    /// Clone the reaction values of another compartment into this one
    void cloneReactions (Compartment *target) const {
        assert(target->numberOfReactions()==0);
        for(auto &r : _internal_reactions.reactions()){
            target->addInternalReactionUnique(unique_ptr<ReactionBase>(
                                            r->clone(target->_species)));
        }
    }

    /// Clone both the species and compartments into this compartment
    void cloneSpeciesReactions(Compartment* C) {
        if(_activated) this->cloneSpecies(C);
        this->cloneReactions(C);
        C->_diffusion_rates = this->_diffusion_rates;
        
    }
    
    /// Clone a compartment
    /// @note - this does not clone the neighbors, just reactions and species
    virtual Compartment* clone() {
        Compartment *C = new Compartment(*this);
        return C;
    }
    
    /// Generate diffusion reactions between this compartment and another
    ///@return - a vector of reactionbases that was just added 
    vector<ReactionBase*> generateDiffusionReactions(Compartment* C);

    /// Generate all diffusion reactions for this compartment and its neighbors
    void generateAllDiffusionReactions();
    
    /// Gives the number of neighbors to this compartment
    size_t numberOfNeighbours() const {return _neighbours.size();}
    
    
    ///Get the species container vector
    SpeciesPtrContainerVector& getSpeciesContainer() {return _species;}
    const SpeciesPtrContainerVector& getSpeciesContainer() const {return _species;}
    
    ///Get the internal reaction container vector
    ReactionPtrContainerVector& getInternalReactionContainer() {return _internal_reactions;}
    const ReactionPtrContainerVector& getInternalReactionContainer() const {return _internal_reactions;}
    
    ///Get the diffusion reaction container vector
    ReactionPtrContainerVector& getDiffusionReactionContainer() {return _diffusion_reactions;}
    const ReactionPtrContainerVector& getDiffusionReactionContainer() const {return _diffusion_reactions;}
    
    /// Get the vector list of neighbors to this compartment
    vector<Compartment*>& getNeighbours() {return _neighbours;}
    const vector<Compartment*>& getNeighbours() const {return _neighbours;}
    
    /// Print the species in this compartment
    void printSpecies() {_species.printSpecies();}
    /// Print the reactions in this compartment
    void printReactions() {
        _internal_reactions.printReactions();
        _diffusion_reactions.printReactions();
    }
    
    /// Check if all species are unique pointers
    bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}
    
    /// Adds the reactions of this compartment to the ChemSim object
    /// @param - chem, a ChemSim object that runs the reaction-diffusion algorithm
    virtual void addChemSimReactions() {
        for(auto &r : _internal_reactions.reactions()) ChemSim::addReaction(r.get());
        for(auto &r : _diffusion_reactions.reactions()) ChemSim::addReaction(r.get());
    }
    
    /// Print properties of this compartment
    virtual void printSelf() {
        cout << this->getFullName() << "\n"
        << "Number of neighbors: " << numberOfNeighbours() << endl;
        printSpecies();
        cout << "Reactions:" << endl;
        printReactions();
    }

};
#endif
