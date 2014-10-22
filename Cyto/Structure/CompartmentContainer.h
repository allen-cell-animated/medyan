//
//  CompartmentContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CompartmentContainer__
#define __CytoSim__CompartmentContainer__

#include <iostream>

#include "common.h"
#include "Compartment.h"

///Key to access instance of CompartmentGrid
class CompartmentGridKey {friend class ChemInitializerImpl;
                          friend class GController;
                          friend class ReactionFilamentTemplate;
#ifdef TESTING
                          public:
#endif //TESTING
                          CompartmentGridKey(){} public: ~CompartmentGridKey(){} };


/// CompartmentGrid class is a simple n-dimensional grid of CompartmentSpatial objects (singleton)

/*!
 *  The CompartmentGrid class is a singleton grid of CompartmentSpatial objects, each of them seperately
 *  holding internal and diffusion reactions, species information, as well as spatial information.
 *  This class is n-dimensional, and the dimension is specified at runtime.
 *
 *  All compartments within CompartmentGrid are indexed by the geometry controller, and this class 
 *  is responsible for assignment of compartment neighbors upon initialization. All initialization
 *  of the CompartmentGrid should be done through GController::initializeGrid().
 *
 *  The _prototype_compartment is used such that all reactions and species can be added to 
 *  the prototype, and this configuration can be copied to all compartments in the grid.
 *  An example of such is below:
 *
 *  @code
 *  GController g;
 *  g.initializeGrid(3, {50, 50, 50}, {5000.0, 5000.0, 5000.0});
 *
 *  Compartment &Cproto = CompartmentGrid::Instance(CompartmentGridKey())->getProtoCompartment();
 *  Species *M1 = Cproto.addSpeciesDiffusing("Myosin",1U);
 *  Cproto.setDiffusionRate(M1,2000);
 *  Species *M2 = Cproto.addSpeciesDiffusing("Fascin",6U);
 *  Cproto.setDiffusionRate(M2,2000);
 *
 *  Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
 *  Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
 *  @endcode
 *
 */
class CompartmentGrid : public Composite {
private:
    Compartment _prototype_compartment; ///< prototype compartment, to be configured before initialization
    SpeciesPtrContainerVector _bulkSpecies; ///<Bulk species in this grid
    ReactionPtrContainerVector _bulkReactions; ///< Bulk reactions in this grid
    
    int _numCompartments; ///<num compartments in the grid
    
    static CompartmentGrid* _instance; ///singleton instance
    
    //private constructor
    CompartmentGrid(int numCompartments) : _numCompartments(numCompartments) {
        
        ///add children
        for(size_t i=0; i<numCompartments; ++i)
        {
            addChild(std::unique_ptr<Component>(new Compartment()));
        }
    }
public:
    
    /// Copying is not allowed
    CompartmentGrid(const CompartmentGrid &rhs) = delete;
    
    /// Assignment is not allowed
    CompartmentGrid& operator=(CompartmentGrid &rhs) = delete;
    
    ///set instance of grid (should only be done at beginning of program)
    ///@note if called, a completely new grid will be created, and the old one erased.
    static void setInstance(CompartmentGridKey k, int numCompartments);
    
    ///Get instance of grid
    static CompartmentGrid* Instance(CompartmentGridKey k);
    
    /// Get name of this compartment grid
    virtual std::string getFullName() const {return std::string("CompartmentGrid");};
    
    ///Activate all compartments
    void activateAll()
    {
        for (auto &c : children())
            static_cast<Compartment*>(c.get())->activate();
        
    }
    
    /// Get the protocompartment from this grid, in order to configure and then initialize
    Compartment& getProtoCompartment() {return _prototype_compartment;}
    const Compartment& getProtoCompartment() const {return _prototype_compartment;}
    
    
    /// Add reactions to all compartments in the grid
    /// @param - chem, a ChemSim object that controls the reaction algorithm
    virtual void addChemSimReactions()
    {
        for(auto &c : children())
        {
            Compartment*C = static_cast<Compartment*>(c.get());
            C->addChemSimReactions();
        }
        
        for(auto &r : _bulkReactions.reactions())
            ChemSim::addReaction(ChemSimReactionKey(), r.get());
        
    }
    
    /// Print properties of this grid
    virtual void printSelf()
    {
        std::cout << getFullName() << std::endl;
        std::cout << "Number of Compartment objects: " << numberOfChildren() << std::endl;
        for(auto &c : children())
            c->printSelf();
    }
    
    ///Add a bulk species to this grid
    template<typename ...Args>
    SpeciesBulk* addSpeciesBulk (Args&& ...args) {
        _bulkSpecies.addSpecies<SpeciesBulk>(std::forward<Args>(args)...);
        return static_cast<SpeciesBulk*>(_bulkSpecies.findSpeciesByIndex(_bulkSpecies.size() - 1));
    }
    
    ///Remove bulk species
    void removeSpeciesBulk(const std::string& name) {_bulkSpecies.removeSpecies(name);}
    
    ///Bulk species finder functions
    SpeciesBulk* findSpeciesBulkByName(const std::string& name) {
        return static_cast<SpeciesBulk*>(_bulkSpecies.findSpeciesByName(name));
    }
    SpeciesBulk* findSpeciesBulkByMolecule(int molecule) {
        return static_cast<SpeciesBulk*>(_bulkSpecies.findSpeciesByMolecule(molecule));
    }
    
    /// Add a bulk reaction to this compartment grid
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addBulkReaction (Args&& ...args)
    {
        //            std::cout << "Compartment::addReaction()..." << std::endl;
        ReactionBase *r = _bulkReactions.addReaction<M,N>(std::forward<Args>(args)...);
        r->setParent(this);
        return r;
    }
    
    /// Add an bulk reaction to this compartment grid
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addBulk(std::initializer_list<Species*> species, float rate)
    {
        ReactionBase *r = _bulkReactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }
    
    /// Add a unique bulk reaction pointer to this compartment
    ReactionBase* addBulkReactionUnique (std::unique_ptr<ReactionBase> &&reaction)
    {
        ReactionBase *r = _bulkReactions.addReactionUnique(std::move(reaction));
        r->setParent(this);
        return r;
    }
    
    /// Remove all bulk reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeBulkReactions (Species* s) {_bulkReactions.removeReactions(s);}

    /// Remove a bulk reaction
    virtual void removeBulkReaction(ReactionBase *r) {_bulkReactions.removeReaction(r);}
    

};


#endif /* defined(__CytoSim__CompartmentContainer__) */
