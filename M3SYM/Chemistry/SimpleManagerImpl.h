
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_SimpleManagerImpl_h
#define M3SYM_SimpleManagerImpl_h

#include <vector>

#include "common.h"

#include "ChemManagerImpl.h"
#include "ReactionManager.h"

//FORWARD DECLARATIONS
class SubSystem;
class InternalFilamentRxnManager;
class CrossFilamentRxnManager;
struct ChemistryData;

/// A concrete implementation of the ChemManagerImpl class.
/// @see ChemManager for documentation on implemented methods.
class SimpleManagerImpl : public ChemManagerImpl {
    
private:
    SubSystem* _subSystem; ///< A pointer to subsytem for creation of callbacks, etc.
    
    /// A list of reactions to add to every new CCylinder
    vector<unique_ptr<InternalFilamentRxnManager>> _IFRxnManagers;
    /// A list of cross filament reactions to add to [CCylinders] (@ref CCylinder)
    vector<unique_ptr<CrossFilamentRxnManager>> _CFRxnManagers;
    
    //@{
    /// Vector of Filament -related species in system
    vector<string> _speciesFilament, _speciesPlusEnd,
                   _speciesMinusEnd, _speciesBound,
                   _speciesLinker, _speciesMotor;
    //@}
    
    //@{
    /// Set up all reaction managers from the setup struct
    void genIFRxnManagers(ChemistryData& chem);
    void genCFRxnManagers(ChemistryData& chem);
    //@}
    
    ///Generate the general, non-filament reactions
    void genGeneralReactions(ChemistryData& chem, Compartment& protoCompartment);
    
    ///Generate bulk reactions
    void genBulkReactions(ChemistryData& chem);
    
    ///Copies species from chem struct
    void copySpecies(ChemistryData& chem);
    
public:
    ///Constructor sets subsystem pointer
    SimpleManagerImpl(SubSystem* subSystem) : _subSystem(subSystem) {}

    virtual void initialize(ChemistryData& chem);
    virtual void initializeCCylinder(CCylinder* cc, Filament* f, bool extensionFront, bool extensionBack, bool creation);
    virtual void updateCCylinder(CCylinder* cc);
    
};

#endif
