
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

#ifndef M3SYM_ChemManager_h
#define M3SYM_ChemManager_h

#include "common.h"

#include "ReactionTemplate.h"
#include "Parser.h"

//FORWARD DECLARATIONS
class Compartment;
class CompartmentGrid;
class ChemSim;

class Cylinder;
class CCylinder;
class CMonomer;

/// For initailizing chemical reactions based on a specific system
/*!
 *  ChemManager is used for initailizing all chemistry in the system.
 *  Initialized by the CController, ChemManager initializes the chemical components
 *  of the CompartmentGrid. The ChemManager can also update chemical components of
 *  a CCylinder using its stored collection of FilamentReactionTemplate, which is
 *  initialized at startup.
 */
class ChemManager {
    
private:
    //@{
    /// Helper functions. Names are pretty self-explanatory.
    void configCMonomer(); void initCMonomer(CMonomer* m, Compartment* c);
    
    void genSpecies(Compartment& protoCompartment);
    
    void genGeneralReactions(Compartment& protoCompartment);
    void genBulkReactions();
    
    void genNucleationReactions();
    void genFilReactionTemplates();
    void genFilBindingManagers();
    //@}
    
public:
    ///Constructor sets subsystem pointer, and loads chem data
    ChemManager(SubSystem* subSystem, ChemistryData chem)
        : _subSystem(subSystem), _chemData(chem) {}
    
    /// Initialize the system, based on the given simulation
    /// Adds all necessary reactions to the ChemSim object
    virtual void initializeSystem(ChemSim* chem);
    
    ///Initializer for chem cylinders, based on the given simulation
    virtual void initializeCCylinder(CCylinder* cc,
                                     bool extensionFront,
                                     bool extensionBack,
                                     bool initialization);
    
    /// Update the copy numbers of all species in the chemical network
    /// @note - this only sets the copy number if the simulation time
    ///         tau has passed the release time of the molecule. This
    ///         function is called at every set of chemical steps to check
    ///         if molecules should be released at the current time.
    virtual void updateCopyNumbers();
    
private:
    //DATA MEMBERS
    SubSystem* _subSystem;   ///< A pointer to subsytem for creation of callbacks, etc.
    ChemistryData _chemData; ///<The chemistry data for the system
    
    /// A list of reactions to add to every new CCylinder
    vector<unique_ptr<FilamentReactionTemplate>> _filRxnTemplates;
};


#endif