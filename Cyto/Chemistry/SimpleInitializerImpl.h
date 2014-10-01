//
//  SimpleInitializerImpl.h
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

///SimpleInitializer is a concrete implementation of the ChemInitailizerImpl class
class SimpleInitializerImpl : public ChemInitializerImpl {

private:
    std::vector<std::unique_ptr<ReactionFilamentTemplate>> _reactionFilamentTemplates; ///< list of reactions to add to every new CCylinder
    
    ///Vectors of all filament, bound, and end species
    std::vector<std::string> _speciesFilament;
    std::vector<std::string> _speciesBound;
    std::vector<std::string> _speciesPlusEnd;
    std::vector<std::string> _speciesMinusEnd;
    
    ///Set up all filament reaction templates from chemsetup struct
    void createFilamentReactionTemplates(ChemistrySpeciesAndReactions& chemSR);
    
public:
    ///initialize the chemical reaction templates and species in this system
    ///@param chemSetup - chemistry setup struct from parsed input file
    virtual void initialize(ChemistrySpeciesAndReactions& chemSR);

    ///Initializer
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    ///@note when initializing, the filaments are filled with the first species listed in the speciesFilament
    /// vector. The active plus and minus end is set to be the first listed as well.
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack);
    
    ///Remove a CCylinder
    virtual void removeCCylinder(Filament* pf, bool retractionFront, bool retractionBack);
    
};

#endif /* defined(__Cyto__SimpleInitializerImpl__) */
