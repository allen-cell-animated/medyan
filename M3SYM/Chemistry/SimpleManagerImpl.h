
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
#include "ReactionManagerImpl.h"
#include "Parser.h"

//FORWARD DECLARATIONS
class SubSystem;
class CMonomer;
class Compartment;

class InternalFilamentRxnManager;
class CrossFilamentRxnManager;

/// A concrete implementation of the ChemManagerImpl class.
/// @see ChemManager for documentation on implemented methods.
class SimpleManagerImpl : public ChemManagerImpl {
    
private:
    ChemistryData _chemData; ///<The chemistry data for the system
    SubSystem* _subSystem; ///< A pointer to subsytem for creation of callbacks, etc.
    
    /// A list of reactions to add to every new CCylinder
    vector<unique_ptr<InternalFilamentRxnManager>> _IFRxnManagers;
    /// A list of cross filament reactions to add to [CCylinders] (@ref CCylinder)
    vector<unique_ptr<CrossFilamentRxnManager>> _CFRxnManagers;
    
    
    /// Configure memory and parameters of CMonomer class
    void configureCMonomer();
    /// Intialize a CMonomer based on system chemistry
    void initCMonomer(CMonomer* m, Compartment* c);
    
    /// Generate the general, non-filament reactions
    void genGeneralReactions(Compartment& protoCompartment);
    /// Generate bulk reactions
    void genBulkReactions();
    /// Generate reactions that create new filaments from
    /// diffusing and/or bulk species
    void genNucleationReactions();
    
    //@{
    /// Set up all Filament ReactionManagers from the setup struct
    void genIFRxnManagers();
    void genCFRxnManagers();
    //@}
    
    /// Generate all species, bulk and diffusing
    void genSpecies(Compartment& protoCompartment);
    
public:
    ///Constructor sets subsystem pointer
    SimpleManagerImpl(SubSystem* subSystem, ChemistryData chem)
        : _subSystem(subSystem), _chemData(chem) {}

    virtual void initialize();
    
    virtual void initializeCCylinder(CCylinder* cc, Filament* f,
                                     bool extensionFront,
                                     bool extensionBack,
                                     bool creation);
    
    virtual void updateCCylinder(CCylinder* cc);
    
    virtual void updateCopyNumbers();
    
};

#endif
