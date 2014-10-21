//
//  InitializerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__SimpleInitializerImpl__
#define __Cyto__SimpleInitializerImpl__

#include <iostream>

#include "ChemInitializerImpl.h"
#include "ReactionTemplate.h"
#include "CMonomer.h"
#include "Parser.h"
#include "common.h"

class SubSystem;

///InitializerImpl is a concrete implementation of the ChemInitailizerImpl class
class InitializerImpl : public ChemInitializerImpl {
    
private:
    SubSystem* _subSystem; ///< ptr to subsytem for creation of callbacks, etc
    
    std::vector<std::unique_ptr<ReactionFilamentTemplate>> _filamentReactionTemplates; ///< list of reactions to add to every new CCylinder
    std::vector<std::unique_ptr<ReactionCrossFilamentTemplate>> _crossFilamentReactionTemplates; ///<list of cross filament reactions to add to CCylinders
    
    ///Vectors of all filament-related species in system
    std::vector<std::string> _speciesFilament, _speciesPlusEnd, _speciesMinusEnd, _speciesBound, _speciesLinker, _speciesMotor;
    
    ///Set up all reaction templates from chemsetup struct
    void generateFilamentReactionTemplates(ChemistrySpeciesAndReactions& chemSR);
    void generateCrossFilamentReactionTemplates(ChemistrySpeciesAndReactions& chemSR);
    
    ///Generate the general, non-filament reactions
    void generateGeneralReactions(ChemistrySpeciesAndReactions& chemSR, Compartment& protoCompartment);
    ///Generate bulk reactions
    void generateBulkReactions(ChemistrySpeciesAndReactions& chemSR);
    
public:
    InitializerImpl(SubSystem* subSystem) : _subSystem(subSystem) {}
    
    ///initialize the chemical reaction templates and species in this system
    ///@param chemSetup - chemistry setup struct from parsed input file
    virtual void initialize(ChemistrySpeciesAndReactions& chemSR);

    ///Initializer
    ///@param length - starting length of the CCylinder initialized
    ///@note when initializing, the filaments are filled with the first species listed in the speciesFilament
    /// vector. The active plus and minus end is set to be the first listed as well.
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack, bool init);
    
    ///add/update cross cylinder reactions that are within range
    virtual void updateCCylinder(CCylinder* cc);
    
};

#endif /* defined(__Cyto__InitializerImpl__) */
