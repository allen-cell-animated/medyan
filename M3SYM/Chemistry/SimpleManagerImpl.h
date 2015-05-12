
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

#ifndef M3SYM_SimpleManagerImpl_h
#define M3SYM_SimpleManagerImpl_h

#include <vector>

#include "common.h"

#include "ChemManagerImpl.h"
#include "ReactionTemplate.h"
#include "Parser.h"

//FORWARD DECLARATIONS
class SubSystem;
class CMonomer;
class Compartment;


/// A concrete implementation of the ChemManagerImpl class.
/// @see ChemManager for documentation on implemented methods.
class SimpleManagerImpl : public ChemManagerImpl {
    
private:
//HELPER FUNCTIONS
    
    /// Configure memory and parameters of CMonomer class
    void configureCMonomer();
    /// Intialize a CMonomer based on system chemistry
    void initCMonomer(CMonomer* m, Compartment* c);
    
    void genGeneralReactions(Compartment& protoCompartment);
    void genBulkReactions();
    
    /// Generate reactions that create new filaments from
    /// diffusing and/or bulk species
    void genNucleationReactions();
    
    /// Set up all [FilamentReactionTemplate](@ref FilamentReactionTemplate)
    /// from the setup struct
    void genFilRxnTemplates();
    /// Set up all [FilamentBindingManagers](@ref FilamentBindingManager)
    /// from the setup struct. Adds to each Compartment.
    void genFilBindingManagers();
    
    void genSpecies(Compartment& protoCompartment);
    
public:
    ///Constructor sets subsystem pointer
    SimpleManagerImpl(SubSystem* subSystem, ChemistryData chem)
    
        : _subSystem(subSystem), _chemData(chem) {}

    virtual void initializeSystem();
    
    virtual void initializeCCylinder(CCylinder* cc, Filament* f,
                                     bool extensionFront,
                                     bool extensionBack,
                                     bool creation);
    
    virtual void updateCopyNumbers();
    
private:
    //DATA MEMBERS
    SubSystem* _subSystem;   ///< A pointer to subsytem for creation of callbacks, etc.
    ChemistryData _chemData; ///<The chemistry data for the system
    
    /// A list of reactions to add to every new CCylinder
    vector<unique_ptr<FilamentReactionTemplate>> _filRxnTemplates;

    
};

#endif
