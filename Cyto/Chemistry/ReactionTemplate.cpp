//
//  ReactionTemplate.cpp
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ReactionTemplate.h"
#include "CompartmentContainer.h"
#include "ChemCallbacks.h"
#include "Bead.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

SubSystem* ReactionFilamentTemplate::_ps = 0;
SubSystem* ReactionCrossFilamentTemplate::_ps = 0;


void PolymerizationPlusEndTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();

    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i+1);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species

        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        ///FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                                 findSpeciesBulkByMolecule(speciesInt));

        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
    
        ///SECOND REACTANT MUST BE PLUS END
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
    
        reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
    
        ///FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);

        productSpecies.push_back(m1->speciesFilament(speciesInt));

    
        ///SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);

        productSpecies.push_back(m2->speciesBound(speciesInt));
    
    
        ///THIRD PRODUCT MUST BE PLUS END
        p = _products[2];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
    
        productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

        ///callbacks
        FilamentExtensionFrontCallback extCallback(pf);
        FilamentPolymerizationFrontCallback polyCallback(pf);
        
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 3>(species, _rate);

        if(i == maxlength - 2)
            boost::signals2::shared_connection_block rcb(rxn->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATION);
    }
}


void PolymerizationMinusEndTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = maxlength - 1; i > 0; i--) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i-1);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        
        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        ///FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                      findSpeciesBulkByMolecule(speciesInt));
        
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        ///SECOND REACTANT MUST BE MINUS END
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        
        ///FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m1->speciesFilament(speciesInt));
        
        
        ///SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesBound(speciesInt));
        
        
        ///THIRD PRODUCT MUST BE MINUS END
        p = _products[2];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
        
        
        FilamentExtensionBackCallback extCallback(pf);
        FilamentPolymerizationBackCallback polyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
        
        if(i == 1)
            boost::signals2::shared_connection_block rcb(rxn->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATION);
    }

}



void PolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species

    
    auto r = _reactants[0];
    SpeciesType type = std::get<1>(r);
    int speciesInt = std::get<0>(r);
    
    ///FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (type == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                  findSpeciesBulkByMolecule(speciesInt));
    
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    ///SECOND REACTANT MUST BE PLUS END
    r = _reactants[1];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
    
    ///FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m1->speciesFilament(speciesInt));
    
    
    ///SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesBound(speciesInt));
    
    
    ///THIRD PRODUCT MUST BE PLUS END
    p = _products[2];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

    
    FilamentPolymerizationFrontCallback polyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATION);
}

void PolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {
    
    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    auto r = _reactants[0];
    SpeciesType type = std::get<1>(r);
    int speciesInt = std::get<0>(r);
    
    ///FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (type == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                  findSpeciesBulkByMolecule(speciesInt));
    
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    ///SECOND REACTANT MUST BE MINUS END
    r = _reactants[1];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
    
    
    ///FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m1->speciesFilament(speciesInt));
    
    
    ///SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesBound(speciesInt));
    
    
    ///THIRD PRODUCT MUST BE MINUS END
    p = _products[2];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));

    FilamentPolymerizationBackCallback polyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATION);
}


void DepolymerizationPlusEndTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();

    ///loop through all monomers
    for(int i = maxlength - 1; i > 0; i--) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i-1);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        
        ///FIRST REACTANT  MUST BE FILAMENT
        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m2->speciesFilament(speciesInt));
        
        ///SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        ///THIRD REACTANT MUST BE PLUSEND
        r = _reactants[2];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);

        reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));

        
        ///FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        if( type == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                      findSpeciesBulkByMolecule(speciesInt));
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        ///SECOND PRODUCT SPECIES MUST BE PLUS END
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesPlusEnd(speciesInt));
        
        FilamentRetractionFrontCallback retCallback(pf);
        FilamentDepolymerizationFrontCallback depolyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        
        if(i == maxlength - 1)
            boost::signals2::shared_connection_block rcb(rxn->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block rcb(rxn->connect(depolyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
    }


}

void DepolymerizationMinusEndTemplate::addReaction(CCylinder* cc, Filament* pf) {

    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i+1);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
 
        ///FIRST REACTANT  MUST BE FILAMENT
        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m2->speciesFilament(speciesInt));
        
        ///SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        ///THIRD REACTANT MUST BE MINUSEND
        r = _reactants[2];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        
        ///FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        if( type == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                     findSpeciesBulkByMolecule(speciesInt));
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        ///SECOND PRODUCT SPECIES MUST BE MINUSEND
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
        
        FilamentRetractionBackCallback retCallback(pf);
        FilamentDepolymerizationBackCallback depolyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        
        if(i == 0)
            boost::signals2::shared_connection_block rcb(rxn->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block rcb(rxn->connect(depolyCallback,false));
            
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
    }
}

void DepolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {
    

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    
    ///FIRST REACTANT  MUST BE FILAMENT
    auto r = _reactants[0];
    SpeciesType type = std::get<1>(r);
    int speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
    
    ///SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesBound(speciesInt));
    
    ///THIRD REACTANT MUST BE PLUSEND
    r = _reactants[2];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
    
    ///FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    if( type == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                 findSpeciesBulkByMolecule(speciesInt));
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    ///SECOND PRODUCT SPECIES MUST BE PLUS END
    p = _products[1];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

    
    FilamentDepolymerizationFrontCallback depolyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(depolyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATION);

}

void DepolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    ///FIRST REACTANT  MUST BE FILAMENT
    auto r = _reactants[0];
    SpeciesType type = std::get<1>(r);
    int speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
    
    ///SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesBound(speciesInt));
    
    ///THIRD REACTANT MUST BE MINUSEND
    r = _reactants[2];
    type = std::get<1>(r);
    speciesInt = std::get<0>(r);
    
    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
    
    
    ///FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    if( type == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                 findSpeciesBulkByMolecule(speciesInt));
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    ///SECOND PRODUCT SPECIES MUST BE MINUSEND
    p = _products[1];
    type = std::get<1>(p);
    speciesInt = std::get<0>(p);
    
    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
    
    FilamentDepolymerizationBackCallback depolyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(depolyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
}


void BasicBindingTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
    
        ///FIRST REACTANT SHOULD BE BOUND
        auto r = _reactants[0];
        auto type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        ///SECOND REACTANT MUST BE BULK OR DIFFUSING
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        if (type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                      findSpeciesBulkByMolecule(speciesInt));
        
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        ///FIRST PRODUCT MUST BE BOUND SPECIES
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m1->speciesBound(speciesInt));
        
    
        ///Add the reaction
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 1>(species, _rate);
        
        //boost::signals2::shared_connection_block rcb(rxn->connect(bindingCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::BASICBINDING);
    }
}


void UnbindingTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        SpeciesType boundType; ReactionType reactionType;
        
        CMonomer* m1 = cc->getCMonomer(i);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        UnbindingCallback unbindingCallback(nullptr, _ps);
        
        ///FIRST REACTANT SHOULD BE LINKER, MOTOR, OR BOUND
        auto r = _reactants[0];
        auto type = std::get<1>(r);
        int speciesInt = std::get<0>(r);

        boundType = type;
        
        if(type == SpeciesType::LINKER) {
            reactantSpecies.push_back(m1->speciesLinker(speciesInt));
            unbindingCallback._s1 = m1->speciesLinker(speciesInt);
            reactionType = ReactionType::LINKERUNBINDING;
        }
        else if(type == SpeciesType::MOTOR) {
            reactantSpecies.push_back(m1->speciesMotor(speciesInt));
            unbindingCallback._s1 = m1->speciesMotor(speciesInt);
            reactionType = ReactionType::MOTORUNBINDING;
        }
        else {
            reactantSpecies.push_back(m1->speciesBound(speciesInt));
            unbindingCallback._s1 = m1->speciesBound(speciesInt);
            reactionType = ReactionType::BASICUNBINDING;
        }
        
        ///FIRST PRODUCT SHOULD BE BULK OR DIFFUSING
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        if(type == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                     findSpeciesBulkByMolecule(speciesInt));
        
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        ///SECOND PRODUCT SHOULD BE BOUND
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m1->speciesBound(speciesInt));
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<1, 2>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(unbindingCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(reactionType);
    }   
}


void LinkerBindingTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    ///check if reaction should be added, i.e. if cc1 and cc2 are within reaction range
    auto c1b1 = cc1->getCylinder()->GetFirstBead()->coordinate;
    auto c1b2 = cc1->getCylinder()->GetSecondBead()->coordinate;
    auto cc1Position = MidPointCoordinate(c1b1, c1b2, 0.5);
    
    auto c2b1 = cc2->getCylinder()->GetFirstBead()->coordinate;
    auto c2b2 = cc2->getCylinder()->GetSecondBead()->coordinate;
    auto cc2Position = MidPointCoordinate(c2b1, c2b2, 0.5);
    
    double dist = TwoPointDistance(cc1Position, cc2Position);
    
    if(_rMin <= dist && dist <= _rMax && cc1->getCylinder()->getFilament() != cc2->getCylinder()->getFilament()) {

        ///Add reaction to middle of both cylinders
        int i, j;
        i = j = int(0.5 * SystemParameters::Geometry().cylinderIntSize);
        
        CMonomer* m1 = cc1->getCMonomer(i);
        CMonomer* m2 = cc2->getCMonomer(j);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species

        ///FIRST AND SECOND REACTANT SHOULD BE BOUND
        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));

        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m2->speciesBound(speciesInt));
        
        ///THIRD REACTANT SHOULD BE BULK OR DIFFUSING
        
        r = _reactants[2];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);

        if(type == SpeciesType::BULK)
           reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                     findSpeciesBulkByMolecule(speciesInt));

        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc1->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }

        
        ///FIRST AND SECOND PRODUCT SHOULD BE LINKER
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        short linkerNumber = speciesInt;
        
        productSpecies.push_back(m1->speciesLinker(speciesInt));
        
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesLinker(speciesInt));

        ///set up callbacks
        LinkerBindingCallback lcallback(cc1->getCylinder(), cc2->getCylinder(), linkerNumber, i, j, _ps);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        rxn->setRMin(_rMin);
        rxn->setRMax(_rMax);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(lcallback,false));
        
        cc1->addCrossCylinderReaction(cc2, rxn);
        rxn->setReactionType(ReactionType::LINKERBINDING);
        
    }
}

void MotorBindingTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    ///check if reaction should be added, i.e. if cc1 and cc2 are within reaction range
    auto c1b1 = cc1->getCylinder()->GetFirstBead()->coordinate;
    auto c1b2 = cc1->getCylinder()->GetSecondBead()->coordinate;
    auto cc1Position = MidPointCoordinate(c1b1, c1b2, 0.5);
    
    auto c2b1 = cc2->getCylinder()->GetFirstBead()->coordinate;
    auto c2b2 = cc2->getCylinder()->GetSecondBead()->coordinate;
    auto cc2Position = MidPointCoordinate(c2b1, c2b2, 0.5);
    
    double dist = TwoPointDistance(cc1Position, cc2Position);
    
    if(_rMin <= dist && dist <= _rMax && cc1->getCylinder()->getFilament() != cc2->getCylinder()->getFilament()) {
        
        ///Add reaction to middle of both cylinders
        int i, j;
        i = j = int(0.5 * SystemParameters::Geometry().cylinderIntSize);
        
        CMonomer* m1 = cc1->getCMonomer(i);
        CMonomer* m2 = cc2->getCMonomer(j);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        
        ///FIRST AND SECOND REACTANT SHOULD BE BOUND
        auto r = _reactants[0];
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        r = _reactants[1];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        reactantSpecies.push_back(m2->speciesBound(speciesInt));
        
        ///THIRD REACTANT SHOULD BE BULK OR DIFFUSING
        
        r = _reactants[2];
        type = std::get<1>(r);
        speciesInt = std::get<0>(r);
        
        if(type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::Instance(compartmentGridKey())->
                                      findSpeciesBulkByMolecule(speciesInt));
        
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc1->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        
        ///FIRST AND SECOND PRODUCT SHOULD BE MOTOR
        auto p = _products[0];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        int motorNumber = speciesInt;
        
        productSpecies.push_back(m1->speciesMotor(speciesInt));
        
        p = _products[1];
        type = std::get<1>(p);
        speciesInt = std::get<0>(p);
        
        productSpecies.push_back(m2->speciesMotor(speciesInt));
        
        ///set up callbacks
        MotorBindingCallback mcallback(cc1->getCylinder(), cc2->getCylinder(), motorNumber, i, j, _ps);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        rxn->setRMin(_rMin);
        rxn->setRMax(_rMax);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(mcallback,false));
        
        cc1->addCrossCylinderReaction(cc2, rxn);
        rxn->setReactionType(ReactionType::MOTORBINDING);


    }
}



