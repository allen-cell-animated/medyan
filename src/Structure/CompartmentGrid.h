
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CompartmentGrid_h
#define MEDYAN_CompartmentGrid_h

#include "common.h"

#include "Compartment.h"

//FORWARD DECLARATIONS
class ChemSim;

/// A simple n-dimensional grid of Compartment objects.

/*!
 *  The CompartmentGrid class is a singleton grid of Compartment objects, each of them 
 *  seperately holding internal and diffusion reactions, species information, as well as
 *  spatial information. This class is n-dimensional, and the dimension is specified at 
 *  runtime.
 *  All compartments within CompartmentGrid are indexed by the GController, and
 *  this class is responsible for assignment of compartment neighbors upon 
 *  initialization. All initialization of the CompartmentGrid should be done through 
 *  GController.initializeGrid().
 *
 *  The _prototype_compartment is used such that all reactions and species can be added 
 *  to the prototype, and this configuration can be copied to all compartments in the grid.
 */
class CompartmentGrid : public Composite {
private:
    Compartment _prototype_compartment; ///< Prototype compartment, to be configured
                                        ///< before initialization
    
    SpeciesPtrContainerVector _bulkSpecies;    ///< Bulk species in this grid
    ReactionPtrContainerVector _bulkReactions; ///< Bulk reactions in this grid
    
public:
    /// Constructor, creates a number of Compartment instances
    CompartmentGrid(int numCompartments) {
        
        //add children
        for(size_t i=0; i<numCompartments; ++i)
            addChild(unique_ptr<Component>(new Compartment()));
    }
    
    /// Get compartments that this grid holds
    vector<Compartment*> getCompartments() {
        
        vector<Compartment*> compartments;
        
        for (auto &c : children())
            compartments.push_back((Compartment*)c.get());
        
        return compartments;
    }
    
    /// Get a compartment at certain index
    Compartment* getCompartment(int index) {
        
        return (Compartment*)(children()[index].get());
    }
    
    /// Set all compartments as active. Used at initialization
    void setAllAsActive() {
        for (auto C : getCompartments()) C->setAsActive();
    }
    
    /// Get name of this compartment grid
    virtual string getFullName() const {return string("CompartmentGrid");};
    
    /// Get the protocompartment from this grid, in order to configure and then initialize
    Compartment& getProtoCompartment() {return _prototype_compartment;}
    const Compartment& getProtoCompartment() const {return _prototype_compartment;}
    
    /// Add reactions to all compartments in the grid
    /// @param - chem, a ChemSim object that controls the reaction algorithm
    virtual void addChemSimReactions(ChemSim* chem);
    
    /// Print properties of this grid
    virtual void printSelf()const {
        cout << getFullName() << endl;
        cout << "Number of Compartment objects: " << numberOfChildren() << endl;
        for(auto &c : children())
            c->printSelf();
    }
    
    /// Add a bulk species to this grid
    template<typename ...Args>
    SpeciesBulk* addSpeciesBulk (Args&& ...args) {
        _bulkSpecies.addSpecies<SpeciesBulk>(forward<Args>(args)...);
        return (SpeciesBulk*)(_bulkSpecies.findSpeciesByIndex(_bulkSpecies.size() - 1));
    }
    
    /// Remove bulk species
    void removeSpeciesBulk(const string& name) {_bulkSpecies.removeSpecies(name);}
    
    /// Bulk species finder functions
    SpeciesBulk* findSpeciesBulkByName(const string& name) {
        return (SpeciesBulk*)(_bulkSpecies.findSpeciesByName(name));
    }
    SpeciesBulk* findSpeciesBulkByMolecule(int molecule) {
        return (SpeciesBulk*)(_bulkSpecies.findSpeciesByMolecule(molecule));
    }
    
    /// Add a bulk reaction to this compartment grid
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addBulkReaction (Args&& ...args)
    {
        ReactionBase *r = _bulkReactions.addReaction<M,N>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    /// Add an bulk reaction to this compartment grid
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addBulk(initializer_list<Species*> species, float rate) {
        ReactionBase *r = _bulkReactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }
    
    /// Add a unique bulk reaction pointer to this compartment
    ReactionBase* addBulkReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _bulkReactions.addReactionUnique(move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Remove all bulk reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeBulkReactions (Species* s) {_bulkReactions.removeReactions(s);}

    /// Remove a bulk reaction
    virtual void removeBulkReaction(ReactionBase *r) {_bulkReactions.removeReaction(r);}
    
    
    /// Count the number of diffusing species with a given name
    species_copy_t countDiffusingSpecies(const string& name);
    /// Count the number of bulk species with a given name
    species_copy_t  countBulkSpecies(const string& name);
    
    ///GetType implementation just returns zero (no CompartmentGrid types yet)
    virtual int getType() {return 0;}
    
    
    ///This function gets all diffusing motor species and creates a virtual position within
    ///the simulation volume (dependent on the compartment). This is used for analyzing motor
    ///gradients when species can either be bound or unbound.
    ///@return - a tuple of the motor ID, the type, and two coordinates defining its end positions
    vector<tuple<int, int, vector<double>, vector<double>>> getDiffusingMotors();
    
};

#endif
