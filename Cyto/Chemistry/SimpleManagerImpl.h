//
//  SimpleManagerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__SimpleManagerImpl__
#define __Cyto__SimpleManagerImpl__

#include <iostream>

#include "common.h"

#include "ChemManagerImpl.h"
#include "ReactionManager.h"

///FORWARD DECLARATIONS
class SubSystem;
struct ChemistryData;

class InternalFilamentRxnManager;
class CrossFilamentRxnManager;

///SimpleManagerImpl is a concrete implementation of the ChemInitailizerImpl class
class SimpleManagerImpl : public ChemManagerImpl {
    
private:
    SubSystem* _subSystem; ///< ptr to subsytem for creation of callbacks, etc
    
    ///< list of reactions to add to every new CCylinder
    vector<unique_ptr<InternalFilamentRxnManager>> _IFRxnManagers;
    ///<list of cross filament reactions to add to CCylinders
    vector<unique_ptr<CrossFilamentRxnManager>> _CFRxnManagers;
    
    ///Vectors of all filament-related species in system
    vector<string> _speciesFilament, _speciesPlusEnd,
                   _speciesMinusEnd, _speciesBound,
                   _speciesLinker, _speciesMotor;
    
    ///Set up all reaction managers from chemsetup struct
    void genIFRxnManagers(ChemistryData& chem);
    
    void genCFRxnManagers(ChemistryData& chem);
    
    ///Generate the general, non-filament reactions
    void genGeneralReactions(ChemistryData& chem, Compartment& protoCompartment);
    
    ///Generate bulk reactions
    void genBulkReactions(ChemistryData& chem);
    
public:
    SimpleManagerImpl(SubSystem* subSystem) : _subSystem(subSystem) {}
    
    /// initialize the chemical reaction templates and species in this system
    /// @param chemSetup - chemistry setup struct from parsed input file
    virtual void initialize(ChemistryData& chem);

    /// SimpleManager
    /// @note when initializing, the filaments are filled with the first species listed in
    /// the speciesFilament vector. The active plus and minus end is set to be the first
    /// listed as well.
    virtual void initializeCCylinder(CCylinder* cc, Filament* f, 
                                     bool extensionFront, bool extensionBack, bool creation);
    
    ///add/update cross cylinder reactions that are within range
    virtual void updateCCylinder(CCylinder* cc);
    
};

#endif /* defined(__Cyto__SimpleManagerImpl__) */
