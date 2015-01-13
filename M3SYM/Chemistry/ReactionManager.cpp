
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
        
        //FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (getType(r) == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(getInt(r)));

        else if(getType(r) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
        }
    
        //SECOND REACTANT MUST BE PLUS END
        r = _reactants[1];
        reactantSpecies.push_back(m1->speciesPlusEnd(getInt(r)));
    
        //FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        productSpecies.push_back(m1->speciesFilament(getInt(p)));
    
        //SECOND PRODUCT MUST BE PLUS END
        p = _products[1];
        productSpecies.push_back(m2->speciesPlusEnd(getInt(p)));
        
        //this reaction also marks an empty bound site
        productSpecies.push_back(m1->speciesBound(0));

        //callback
        FilamentPolymerizationFrontCallback polyCallback(cc->getCylinder());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS>(species, _rate);

        boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATIONPLUSEND);
    }
    
    //add extension callback reaction
    CMonomer* m = cc->getCMonomer(maxlength - 1);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    auto r = _reactants[0];
    //FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (getType(r) == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::instance()->
                                  findSpeciesBulkByMolecule(getInt(r)));
    
    else if(getType(r) == SpeciesType::DIFFUSING) {
        Compartment* c = cc->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
    }
    
    //SECOND REACTANT MUST BE PLUS END
    r = _reactants[1];
    reactantSpecies.push_back(m->speciesPlusEnd(getInt(r)));
    
    //FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    productSpecies.push_back(m->speciesFilament(getInt(p)));
    
    //this reaction also marks an empty bound site
    productSpecies.push_back(m->speciesBound(0));
    
    short plusEndProduct = getInt(_products[1]);
    
    //callbacks
    FilamentExtensionFrontCallback extCallback(cc->getCylinder(), plusEndProduct);
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS - 1>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(extCallback,false));
    
    cc->addInternalReaction(rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATIONPLUSEND);
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
        //FIRST REACTANT MUST BE BULK OR DIFFUSING
        if (getType(r) == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(getInt(r)));
        
        else if(getType(r) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
        }
        
        //SECOND REACTANT MUST BE MINUS END
        r = _reactants[1];
        reactantSpecies.push_back(m1->speciesMinusEnd(getInt(r)));
        
        //FIRST PRODUCT MUST BE FILAMENT
        auto p = _products[0];
        productSpecies.push_back(m1->speciesFilament(getInt(p)));
        
        //SECOND PRODUCT MUST BE MINUS END
        p = _products[1];
        productSpecies.push_back(m2->speciesMinusEnd(getInt(p)));
        
        //this reaction also marks an empty bound site
        productSpecies.push_back(m1->speciesBound(0));
        
        FilamentPolymerizationBackCallback polyCallback(cc->getCylinder());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(polyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATIONMINUSEND);
    }
    
    //add the extension callback
    CMonomer* m = cc->getCMonomer(0);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    auto r = _reactants[0];
    
    //FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (getType(r) == SpeciesType::BULK)
        reactantSpecies.push_back(CompartmentGrid::instance()->
                                  findSpeciesBulkByMolecule(getInt(r)));
    
    else if(getType(r)  == SpeciesType::DIFFUSING) {
        Compartment* c = cc->getCompartment();
        reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
    }
    
    //SECOND REACTANT MUST BE MINUS END
    r = _reactants[1];
    reactantSpecies.push_back(m->speciesMinusEnd(getInt(r)));
    
    //FIRST PRODUCT MUST BE FILAMENT
    auto p = _products[0];
    productSpecies.push_back(m->speciesFilament(getInt(p)));
    
    //this reaction also marks an empty bound site
    productSpecies.push_back(m->speciesBound(0));
    
    auto minusEndType = get<0>(_products[1]);
    
    FilamentExtensionBackCallback extCallback(cc->getCylinder(), minusEndType);
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS - 1>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(extCallback,false));
    
    cc->addInternalReaction(rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATIONMINUSEND);
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
        reactantSpecies.push_back(m2->speciesFilament(getInt(r)));
        
        //SECOND REACTANT MUST BE PLUSEND
        r = _reactants[1];
        reactantSpecies.push_back(m1->speciesPlusEnd(getInt(r)));
        
        //this reaction also needs an empty bound site
        reactantSpecies.push_back(m2->speciesBound(0));

        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        if( getType(p) == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                      findSpeciesBulkByMolecule(getInt(p)));
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS END
        p = _products[1];
        productSpecies.push_back(m2->speciesPlusEnd(getInt(p)));
        
        FilamentDepolymerizationFrontCallback depolyCallback(cc->getCylinder());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
            new Reaction<DEPOLYREACTANTS,DEPOLYPRODUCTS>(species, _rate);
        
        boost::signals2::shared_connection_block
            rcb(rxn->connect(depolyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATIONPLUSEND);
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
        reactantSpecies.push_back(m2->speciesFilament(getInt(r)));
        
        //SECOND REACTANT MUST BE MINUSEND
        r = _reactants[1];
        reactantSpecies.push_back(m1->speciesMinusEnd(getInt(r)));
        
        //this reaction also needs an empty bound site
        reactantSpecies.push_back(m2->speciesBound(0));
        
        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        if(getType(p) == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(getInt(p)));
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        //SECOND PRODUCT SPECIES MUST BE MINUSEND
        p = _products[1];
        productSpecies.push_back(m2->speciesMinusEnd(getInt(p)));
        
        FilamentDepolymerizationBackCallback depolyCallback(cc->getCylinder());
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
            new Reaction<DEPOLYREACTANTS,DEPOLYPRODUCTS>(species, _rate);
        
        boost::signals2::shared_connection_block
            rcb(rxn->connect(depolyCallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATIONMINUSEND);
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
    reactantSpecies.push_back(m2->speciesFilament(getInt(r)));
    
    //SECOND REACTANT MUST BE PLUSEND
    r = _reactants[1];
    reactantSpecies.push_back(m1->speciesPlusEnd(getInt(r)));
    
    //this reaction also needs an empty bound site
    reactantSpecies.push_back(m2->speciesBound(0));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(getInt(p)));
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    //SECOND PRODUCT SPECIES MUST BE PLUS END
    p = _products[1];
    productSpecies.push_back(m2->speciesPlusEnd(getInt(p)));

    FilamentRetractionFrontCallback retCallback(cc1->getCylinder());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<DEPOLYREACTANTS,DEPOLYPRODUCTS>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(retCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATIONPLUSEND);
}

void DepolyMinusEndManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    //FIRST REACTANT  MUST BE FILAMENT
    auto r = _reactants[0];
    reactantSpecies.push_back(m2->speciesFilament(getInt(r)));
    
    //SECOND REACTANT MUST BE MINUSEND
    r = _reactants[1];
    reactantSpecies.push_back(m1->speciesMinusEnd(getInt(r)));
    
    //this reaction also needs an empty bound site
    reactantSpecies.push_back(m2->speciesBound(0));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(getInt(p)));
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    //SECOND PRODUCT SPECIES MUST BE MINUSEND
    p = _products[1];
    productSpecies.push_back(m2->speciesMinusEnd(getInt(p)));
    
    FilamentRetractionBackCallback retCallback(cc1->getCylinder());
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn = new Reaction<DEPOLYREACTANTS,DEPOLYPRODUCTS>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(retCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATIONMINUSEND);
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
        int motorType = getInt(r);
        
        //FIRST REACTANT MUST BE MOTOR
        reactantSpecies.push_back(m1->speciesMotor(motorType));
        
        //SECOND REACTANT MUST BE BOUND
        r = _reactants[1];
        int boundType = getInt(r);
        
        reactantSpecies.push_back(m2->speciesBound(boundType));
        
        //FIRST PRODUCT MUST BE MOTOR
        auto p = _products[0];
        productSpecies.push_back(m2->speciesMotor(getInt(p)));
        
        //SECOND PRODUCT MUST BE BOUND
        p = _products[1];
        productSpecies.push_back(m1->speciesBound(getInt(p)));
        
        //callbacks
        MotorWalkingForwardCallback
            motorMoveCallback(cc->getCylinder(), site1, site2,
                              motorType, boundType, _ps);
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
            new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
        
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
    int motorType = getInt(r);
    
    //FIRST REACTANT MUST BE MOTOR
    reactantSpecies.push_back(m1->speciesMotor(motorType));
    
    //SECOND REACTANT MUST BE BOUND
    r = _reactants[1];
    int boundType = getInt(r);
    
    reactantSpecies.push_back(m2->speciesBound(boundType));
    
    //FIRST PRODUCT MUST BE MOTOR
    auto p = _products[0];
    productSpecies.push_back(m2->speciesMotor(getInt(p)));
    
    //SECOND PRODUCT MUST BE BOUND
    p = _products[1];
    productSpecies.push_back(m1->speciesBound(getInt(p)));
    
    //callbacks
    MotorMovingCylinderForwardCallback
        motorChangeCallback(cc1->getCylinder(), cc2->getCylinder(),
                            _bindingSites.back(), motorType, boundType, _ps);
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn =
        new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
    
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
        SpeciesType type = getType(r);
        int speciesInt = getInt(r);
        
        if(type == SpeciesType::FILAMENT)
            reactantSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)
            reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND)
            reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //FIRST PRODUCT MUST BE FILAMENT, PLUS OR MINUS SPECIES
        auto p = _products[0];
        type = getType(p);
        speciesInt = getInt(p);
        
        if(type == SpeciesType::FILAMENT)
            productSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)
            productSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND)
            productSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
            new Reaction<AGINGREACTANTS,AGINGPRODUCTS>(species, _rate);
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::AGING);
    }
}


void DestructionManager::addReaction(CCylinder* cc) {

    //loop through all monomers of filament
    int maxlength = cc->getSize();
    
    //loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        CMonomer* m2 = cc->getCMonomer(i+1);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        //FIRST REACTANT MUST BE PLUS END
        auto r = _reactants[0];
        reactantSpecies.push_back(m2->speciesPlusEnd(getInt(r)));
        
        //SECOND REACTANT MUST BE MINUS END
        r = _reactants[1];
        reactantSpecies.push_back(m1->speciesMinusEnd(getInt(r)));
        
        //ALL PRODUCTS MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        if(getType(p) == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(getInt(p)));
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        p = _products[1];
        
        if(getType(p) == SpeciesType::BULK)
            productSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(getInt(p)));
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        FilamentDestructionCallback dcallback(cc->getCylinder());
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
            new Reaction<DESTRUCTIONREACTANTS,DESTRUCTIONPRODUCTS>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(dcallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::FILAMENTDESTRUCTION);
    }
}

void DestructionManager::addReaction(CCylinder* cc1, CCylinder* cc2) {

    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //FIRST REACTANT MUST BE PLUS END
    auto r = _reactants[0];
    reactantSpecies.push_back(m2->speciesPlusEnd(getInt(r)));
    
    //SECOND REACTANT MUST BE MINUS END
    r = _reactants[1];
    reactantSpecies.push_back(m1->speciesMinusEnd(getInt(r)));
    
    //ALL PRODUCTS MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(getInt(p)));
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    p = _products[1];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(CompartmentGrid::instance()->
                                 findSpeciesBulkByMolecule(getInt(p)));
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    FilamentDestructionCallback dcallback(cc1->getCylinder());
    
    //Add the reaction
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn =
        new Reaction<DESTRUCTIONREACTANTS,DESTRUCTIONPRODUCTS>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(rxn->connect(dcallback,false));
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::FILAMENTDESTRUCTION);
}

void SeveringManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = _bindingSites.begin(); it != _bindingSites.end(); it++) {
        
        int site = *(it);
        CMonomer* m = cc->getCMonomer(site);
        vector<Species*> reactantSpecies;
        
        //REACTANT MUST BE FILAMENT
        auto r = _reactants[0];
        reactantSpecies.push_back(m->speciesFilament(getInt(r)));
        
        //IMPLICITLY NEEDS AN EMPTY BOUND
        reactantSpecies.push_back(m->speciesBound(0));
        
        FilamentSeveringCallback scallback(cc->getCylinder());
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        ReactionBase* rxn =
        new Reaction<SEVERINGREACTANTS + 1,SEVERINGPRODUCTS>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(scallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::SEVERING);
    }
}

void BranchingManager::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = _bindingSites.begin(); it != _bindingSites.end(); it++) {
        
        int site = *(it);
        CMonomer* m = cc->getCMonomer(site);
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        
        //FIRST TWO REACTANTS MUST BE BULK OR DIFFUSING
        auto r = _reactants[0];
        if(getType(r) == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(getInt(r)));
        else if(getType(r) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
        }
        
        r = _reactants[1];
        if(getType(r) == SpeciesType::BULK)
            reactantSpecies.push_back(CompartmentGrid::instance()->
                                     findSpeciesBulkByMolecule(getInt(r)));
        else if(getType(r) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
        }
        
        //THIRD REACTANT MUST BE BOUND
        r = _reactants[2];
        reactantSpecies.push_back(m->speciesBound(getInt(r)));
        
        //FIRST PRODUCT MUST BE BRANCH
        auto p = _products[0];
        int branchType = getInt(p);
        productSpecies.push_back(m->speciesBrancher(branchType));
                                 
        short plusEndProduct = getInt(_products[1]);
        
        BranchingPointCreationCallback
        bcallback(cc->getCylinder(), branchType, plusEndProduct, site, _offRate, _ps);
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
        new Reaction<BRANCHINGREACTANTS,BRANCHINGPRODUCTS - 1>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(rxn->connect(bcallback,false));
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::BRANCHING);
    }
}


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
                reactantSpecies.push_back(m1->speciesBound(getInt(r)));

                r = _reactants[1];
                reactantSpecies.push_back(m2->speciesBound(getInt(r)));
                
                //THIRD REACTANT SHOULD BE BULK OR DIFFUSING
                r = _reactants[2];
                if(getType(r) == SpeciesType::BULK)
                   reactantSpecies.push_back(CompartmentGrid::instance()->
                                             findSpeciesBulkByMolecule(getType(r)));

                else if(getType(r) == SpeciesType::DIFFUSING) {
                    Compartment* c = cc1->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(getType(r)));
                }

                //FIRST AND SECOND PRODUCT SHOULD BE LINKER
                auto p = _products[0];
                short linkerNumber = getInt(p);
                
                productSpecies.push_back(m1->speciesLinker(linkerNumber));
                
                p = _products[1];
                productSpecies.push_back(m2->speciesLinker(getInt(p)));

                //set up callbacks
                LinkerBindingCallback
                    lcallback(cc1->getCylinder(), cc2->getCylinder(),
                              linkerNumber, i, j, _offRate, _ps);
                
                //Add the reaction. If it needs a callback then attach
                vector<Species*> onSpecies = reactantSpecies;
                onSpecies.insert(onSpecies.end(), productSpecies.begin(),
                                 productSpecies.end());
                ReactionBase* onRxn =
                    new Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>
                    (onSpecies, _onRate);
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
                reactantSpecies.push_back(m1->speciesBound(getInt(r)));
                
                r = _reactants[1];
                reactantSpecies.push_back(m2->speciesBound(getInt(r)));
                
                //THIRD REACTANT SHOULD BE BULK OR DIFFUSING
                r = _reactants[2];
                if(getType(r) == SpeciesType::BULK)
                    reactantSpecies.push_back(CompartmentGrid::instance()->
                                              findSpeciesBulkByMolecule(getInt(r)));
                
                else if(getType(r) == SpeciesType::DIFFUSING) {
                    Compartment* c = cc1->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(getInt(r)));
                }
                
                //FIRST AND SECOND PRODUCT SHOULD BE MOTOR
                auto p = _products[0];
                short motorNumber = getInt(p);
                productSpecies.push_back(m1->speciesMotor(motorNumber));
                
                p = _products[1];
                productSpecies.push_back(m2->speciesMotor(getInt(p)));

                //set up callbacks
                MotorBindingCallback
                    mcallback(cc1->getCylinder(), cc2->getCylinder(),
                              motorNumber, i, j, _offRate, _ps);
                
                //Add the reaction. If it needs a callback then attach
                vector<Species*> onSpecies = reactantSpecies;
                onSpecies.insert(onSpecies.end(), productSpecies.begin(),
                                 productSpecies.end());
                ReactionBase* onRxn =
                    new Reaction<LMBINDINGREACTANTS, LMBINDINGPRODUCTS>
                    (onSpecies, _onRate);
                onRxn->setReactionType(ReactionType::MOTORBINDING);
                
                boost::signals2::shared_connection_block
                    rcb(onRxn->connect(mcallback,false));
                
                cc1->addCrossCylinderReaction(cc2, onRxn);
            }
        }
    }
}


