
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

#include "ReactionManager.h"

#include "ChemCallbacks.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

SubSystem* InternalFilamentRxnManager::_ps = 0;
SubSystem* CrossFilamentRxnManager::_ps = 0;

void PolyPlusEndManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers of filament
    int maxlength = cc->getSize();
    
    //loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i+1);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //loop through reactants, products. find all species

        auto r = _reactants[0];
        SpeciesType type = get<1>(r);
        int speciesInt = get<0>(r);
        
        //FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(speciesInt));

        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
    
        //SECOND REACTANT MUST BE PLUS END
        r = _reactants[1];
        speciesInt = get<0>(r);
    
        reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
        //FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        speciesInt = get<0>(p);

        productSpecies.push_back(m1->speciesFilament(speciesInt));

        //SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        speciesInt = get<0>(p);

        productSpecies.push_back(m2->speciesBound(speciesInt));
    
        //THIRD PRODUCT MUST BE PLUS END
        p = _products[2];
        speciesInt = get<0>(p);
    
        productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

        //callbacks
        FilamentExtensionFrontCallback
            extCallback(cc->getCylinder()->getFilament());
        FilamentPolymerizationFrontCallback
            polyCallback(cc->getCylinder()->getFilament());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 3>(species, _rate);

        if(i == maxlength - 2)
            boost::signals2::shared_connection_block
                rcb(rxn->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block
                rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATION);
    }
}


void PolyMinusEndManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers of filament
    int maxlength = cc->getSize();
    
    //loop through all monomers
    for(int i = maxlength - 1; i > 0; i--) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i-1);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //loop through reactants, products. find all species
        
        auto r = _reactants[0];
        SpeciesType type = get<1>(r);
        int speciesInt = get<0>(r);
        
        //FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (type == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(speciesInt));
        
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        //SECOND REACTANT MUST BE MINUS END
        r = _reactants[1];
        speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m1->speciesFilament(speciesInt));
        
        //SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m2->speciesBound(speciesInt));
        
        //THIRD PRODUCT MUST BE MINUS END
        p = _products[2];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
        
        
        FilamentExtensionBackCallback
            extCallback(cc->getCylinder()->getFilament());
        FilamentPolymerizationBackCallback
            polyCallback(cc->getCylinder()->getFilament());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
        
        if(i == 1)
            boost::signals2::shared_connection_block
                rcb(rxn->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block
                rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATION);
    }

}

void PolyPlusEndManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    auto r = _reactants[0];
    SpeciesType type = get<1>(r);
    int speciesInt = get<0>(r);
    
    //FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (type == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::instance()->
                                  findSpeciesBulkByMolecule(speciesInt));
    
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    //SECOND REACTANT MUST BE PLUS END
    r = _reactants[1];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
    //FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m1->speciesFilament(speciesInt));
    
    //SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesBound(speciesInt));

    //THIRD PRODUCT MUST BE PLUS END
    p = _products[2];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

    FilamentPolymerizationFrontCallback
        polyCallback(cc1->getCylinder()->getFilament());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block
        rcb(rxn->connect(polyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATION);
}

void PolyMinusEndManager::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    CMonomer* m1 = cc2->getCMonomer(0);
    CMonomer* m2 = cc1->getCMonomer(cc1->getSize() - 1);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    auto r = _reactants[0];
    SpeciesType type = get<1>(r);
    int speciesInt = get<0>(r);
    
    //FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (type == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::instance()->
                                  findSpeciesBulkByMolecule(speciesInt));
    
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    //SECOND REACTANT MUST BE MINUS END
    r = _reactants[1];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
    
    //FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m1->speciesFilament(speciesInt));
    
    //SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesBound(speciesInt));
    
    //THIRD PRODUCT MUST BE MINUS END
    p = _products[2];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));

    FilamentPolymerizationBackCallback
        polyCallback(cc1->getCylinder()->getFilament());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block
        rcb(rxn->connect(polyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATION);
}


void DepolyPlusEndManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers of filament
    int maxlength = cc->getSize();

    //loop through all monomers
    for(int i = maxlength - 1; i > 0; i--) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i-1);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //loop through reactants, products. find all species
        
        //FIRST REACTANT  MUST BE FILAMENT
        auto r = _reactants[0];
        int speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m2->speciesFilament(speciesInt));
        
        //SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        //THIRD REACTANT MUST BE PLUSEND
        r = _reactants[2];
        speciesInt = get<0>(r);

        reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));

        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        SpeciesType type = get<1>(p);
        speciesInt = get<0>(p);
        
        if( type == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(speciesInt));
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS END
        p = _products[1];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m2->speciesPlusEnd(speciesInt));
        
        FilamentRetractionFrontCallback
            retCallback(cc->getCylinder()->getFilament());
        FilamentDepolymerizationFrontCallback
            depolyCallback(cc->getCylinder()->getFilament());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        
        if(i == maxlength - 1)
            boost::signals2::shared_connection_block
                rcb(rxn->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block
                rcb(rxn->connect(depolyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
    }
}

void DepolyMinusEndManager::addReaction(CCylinder* cc) {

    //loop through all monomers of filament
    int maxlength = cc->getSize();
    
    //loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i+1);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //loop through reactants, products. find all species
 
        //FIRST REACTANT  MUST BE FILAMENT
        auto r = _reactants[0];
        int speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m2->speciesFilament(speciesInt));
        
        //SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m1->speciesBound(speciesInt));
        
        //THIRD REACTANT MUST BE MINUSEND
        r = _reactants[2];
        speciesInt = get<0>(r);
        
        reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        SpeciesType type = get<1>(p);
        speciesInt = get<0>(p);
        
        if( type == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(speciesInt));
        else if(type == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
        }
        
        //SECOND PRODUCT SPECIES MUST BE MINUSEND
        p = _products[1];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
        
        FilamentRetractionBackCallback
            retCallback(cc->getCylinder()->getFilament());
        FilamentDepolymerizationBackCallback
            depolyCallback(cc->getCylinder()->getFilament());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
        
        if(i == 0)
            boost::signals2::shared_connection_block
                rcb(rxn->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block
                rcb(rxn->connect(depolyCallback,false));
            
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
    }
}

void DepolyPlusEndManager::addReaction(CCylinder* cc1, CCylinder* cc2) {
    

    CMonomer* m1 = cc2->getCMonomer(0);
    CMonomer* m2 = cc1->getCMonomer(cc1->getSize() - 1);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    
    //FIRST REACTANT  MUST BE FILAMENT
    auto r = _reactants[0];
    int speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
    
    //SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesBound(speciesInt));
    
    //THIRD REACTANT MUST BE PLUSEND
    r = _reactants[2];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    SpeciesType type = get<1>(p);
    speciesInt = get<0>(p);
    
    if( type == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(speciesInt));
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    //SECOND PRODUCT SPECIES MUST BE PLUS END
    p = _products[1];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));

    FilamentDepolymerizationFrontCallback
        depolyCallback(cc1->getCylinder()->getFilament());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block
        rcb(rxn->connect(depolyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATION);

}

void DepolyMinusEndManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    //FIRST REACTANT  MUST BE FILAMENT
    auto r = _reactants[0];
    int speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
    
    //SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesBound(speciesInt));
    
    //THIRD REACTANT MUST BE MINUSEND
    r = _reactants[2];
    speciesInt = get<0>(r);
    
    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    SpeciesType type = get<1>(p);
    speciesInt = get<0>(p);
    
    if( type == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(speciesInt));
    else if(type == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
    }
    
    //SECOND PRODUCT SPECIES MUST BE MINUSEND
    p = _products[1];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
    
    FilamentDepolymerizationBackCallback
        depolyCallback(cc1->getCylinder()->getFilament());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block
        rcb(rxn->connect(depolyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATION);
}

void MotorWalkFManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = _bindingSites.begin(); it != _bindingSites.end() - 1; it++) {

        int site1 = *(it);
        int site2 = *(it+1);
        
        CMonomer* m1 = cc->getCMonomer(site1);
        CMonomer* m2 = cc->getCMonomer(site2);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //loop through reactants, products. find all species
        auto r = _reactants[0];
        int speciesInt = get<0>(r);
        int motorType = speciesInt;
        
        //FIRST REACTANT MUST BE MOTOR
        reactantSpecies.push_back(m1->speciesMotor(speciesInt));
        
        //SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        speciesInt = get<0>(r);
        int boundType = speciesInt;
        
        reactantSpecies.push_back(m2->speciesBound(speciesInt));
        
        //FIRST PRODUCT MUST BE MOTOR
        auto p = _products[0];
        speciesInt = get<0>(p);
    
        productSpecies.push_back(m2->speciesMotor(speciesInt));
        
        //SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        speciesInt = get<0>(p);
        
        productSpecies.push_back(m1->speciesBound(speciesInt));
        
        //callbacks
        MotorWalkingForwardCallback
            motorMoveCallback(cc->getCylinder(), site1, site2,
                              motorType, boundType, _ps);
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<2, 2>(species, _rate);
        
        boost::signals2::shared_connection_block
            rcb(rxn->connect(motorMoveCallback, false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::MOTORWALKINGFORWARD);
    }
}

void MotorWalkFManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    CMonomer* m1 = cc1->getCMonomer(_bindingSites.back());
    CMonomer* m2 = cc2->getCMonomer(_bindingSites.front());
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;

    //loop through reactants, products. find all species
    auto r = _reactants[0];
    int speciesInt = get<0>(r);
    int motorType = speciesInt;
    
    //FIRST REACTANT MUST BE MOTOR
    reactantSpecies.push_back(m1->speciesMotor(speciesInt));
    
    //SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    speciesInt = get<0>(r);
    int boundType = speciesInt;
    
    reactantSpecies.push_back(m2->speciesBound(speciesInt));
    
    //FIRST PRODUCT MUST BE MOTOR
    auto p = _products[0];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m2->speciesMotor(speciesInt));
    
    //SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    speciesInt = get<0>(p);
    
    productSpecies.push_back(m1->speciesBound(speciesInt));
    
    //callbacks
    MotorMovingCylinderForwardCallback
        motorChangeCallback(cc1->getCylinder(), cc2->getCylinder(),
                            _bindingSites.back(), motorType, boundType, _ps);
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<2, 2>(species, _rate);
    
    boost::signals2::shared_connection_block
        rcb(rxn->connect(motorChangeCallback, false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::MOTORWALKINGFORWARD);
}

void AgingManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers of filament
    int maxlength = cc->getSize();
    
    //loop through all monomers
    for(int i = 0; i <= maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //FIRST REACTANT SHOULD BE FILAMENT, PLUS OR MINUS SPECIES
        auto r = _reactants[0];
        SpeciesType type = get<1>(r);
        int speciesInt = get<0>(r);
        
        if(type == SpeciesType::FILAMENT)
            reactantSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)
            reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND)
            reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //FIRST PRODUCT MUST BE FILAMENT, PLUS OR MINUS SPECIES
        auto p = _products[0];
        type = get<1>(p);
        speciesInt = get<0>(p);
        
        if(type == SpeciesType::FILAMENT)
            productSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)
            productSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND)
            productSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<1, 1>(species, _rate);
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::AGING);
    }
}

void AgingManager::addReaction(CCylinder* cc1, CCylinder* cc2) {}

void LinkerRxnManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    //Add reaction to the binding sites that are within range
    for(int i : _bindingSites) {
        for(int j : _bindingSites) {
            
            //check if these sites are within range
            auto coord1 =
                midPointCoordinate(cc1->getCylinder()->getFirstBead()->coordinate,
                                   cc1->getCylinder()->getSecondBead()->coordinate,
                         (double)i / SystemParameters::Geometry().cylinderIntSize);
            
            auto coord2 =
                midPointCoordinate(cc2->getCylinder()->getFirstBead()->coordinate,
                                   cc2->getCylinder()->getSecondBead()->coordinate,
                         (double)j / SystemParameters::Geometry().cylinderIntSize);
            
            double dist = twoPointDistance(coord1, coord2);
            
            if(dist < _rMax && dist > _rMin) {
            
                CMonomer* m1 = cc1->getCMonomer(i);
                CMonomer* m2 = cc2->getCMonomer(j);
                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;
                
                //loop through reactants, products. find all species

                //FIRST AND SECOND REACTANT SHOULD BE BOUND
                auto r = _reactants[0];
                int speciesInt = get<0>(r);
                
                reactantSpecies.push_back(m1->speciesBound(speciesInt));

                r = _reactants[1];
                speciesInt = get<0>(r);
                
                reactantSpecies.push_back(m2->speciesBound(speciesInt));
                
                //THIRD REACTANT SHOULD BE BULK OR DIFFUSING
                r = _reactants[2];
                SpeciesType type = get<1>(r);
                speciesInt = get<0>(r);

                if(type == SpeciesType::BULK)
                   reactantSpecies.push_back(CompartmentGrid::instance()->
                                             findSpeciesBulkByMolecule(speciesInt));

                else if(type == SpeciesType::DIFFUSING) {
                    Compartment* c = cc1->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                }

                //FIRST AND SECOND PRODUCT SHOULD BE LINKER
                auto p = _products[0];
                speciesInt = get<0>(p);
                short linkerNumber = speciesInt;
                
                productSpecies.push_back(m1->speciesLinker(speciesInt));
                
                p = _products[1];
                speciesInt = get<0>(p);
                
                productSpecies.push_back(m2->speciesLinker(speciesInt));

                //set up callbacks
                LinkerBindingCallback
                    lcallback(cc1->getCylinder(), cc2->getCylinder(),
                              linkerNumber, i, j, _offRate, _ps);
                
                //Add the reaction. If it needs a callback then attach
                vector<Species*> onSpecies = reactantSpecies;
                onSpecies.insert(onSpecies.end(), productSpecies.begin(),
                                 productSpecies.end());
                ReactionBase* onRxn = new Reaction<3, 2>(onSpecies, _onRate);
                onRxn->setReactionType(ReactionType::LINKERBINDING);
                
                boost::signals2::shared_connection_block
                    rcb(onRxn->connect(lcallback,false));
                
                cc1->addCrossCylinderReaction(cc2, onRxn);
            }
        }
    }
}

void MotorRxnManager::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    //Add reaction to the binding sites that are within range
    for(int i : _bindingSites) {
        for(int j : _bindingSites) {
            
            //check if these sites are within range
            auto coord1 =
                midPointCoordinate(cc1->getCylinder()->getFirstBead()->coordinate,
                                   cc1->getCylinder()->getSecondBead()->coordinate,
                         (double)i / SystemParameters::Geometry().cylinderIntSize);
            
            auto coord2 =
                midPointCoordinate(cc2->getCylinder()->getFirstBead()->coordinate,
                                   cc2->getCylinder()->getSecondBead()->coordinate,
                         (double)j / SystemParameters::Geometry().cylinderIntSize);
            
            double dist = twoPointDistance(coord1, coord2);
            
            if(dist < _rMax && dist > _rMin) {
                
                CMonomer* m1 = cc1->getCMonomer(i);
                CMonomer* m2 = cc2->getCMonomer(j);
                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;
                
                //loop through reactants, products. find all species
                
                //FIRST AND SECOND REACTANT SHOULD BE BOUND
                auto r = _reactants[0];
                int speciesInt = get<0>(r);
                
                reactantSpecies.push_back(m1->speciesBound(speciesInt));
                
                r = _reactants[1];
                speciesInt = get<0>(r);
                
                reactantSpecies.push_back(m2->speciesBound(speciesInt));
                
                //THIRD REACTANT SHOULD BE BULK OR DIFFUSING
                r = _reactants[2];
                SpeciesType type = get<1>(r);
                speciesInt = get<0>(r);
                
                if(type == SpeciesType::BULK)
                    reactantSpecies.push_back(CompartmentGrid::instance()->
                                              findSpeciesBulkByMolecule(speciesInt));
                
                else if(type == SpeciesType::DIFFUSING) {
                    Compartment* c = cc1->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                }
                
                //FIRST AND SECOND PRODUCT SHOULD BE MOTOR
                auto p = _products[0];
                speciesInt = get<0>(p);
                short motorNumber = speciesInt;
                
                productSpecies.push_back(m1->speciesMotor(speciesInt));
                
                p = _products[1];
                speciesInt = get<0>(p);
                
                productSpecies.push_back(m2->speciesMotor(speciesInt));

                //set up callbacks
                MotorBindingCallback
                    mcallback(cc1->getCylinder(), cc2->getCylinder(),
                              motorNumber, i, j, _offRate, _ps);
                
                //Add the reaction. If it needs a callback then attach
                vector<Species*> onSpecies = reactantSpecies;
                onSpecies.insert(onSpecies.end(), productSpecies.begin(),
                                 productSpecies.end());
                ReactionBase* onRxn = new Reaction<3, 2>(onSpecies, _onRate);
                onRxn->setReactionType(ReactionType::MOTORBINDING);
                
                boost::signals2::shared_connection_block
                    rcb(onRxn->connect(mcallback,false));
                
                cc1->addCrossCylinderReaction(cc2, onRxn);
            }
        }
    }
}


