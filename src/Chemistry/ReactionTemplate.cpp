
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ReactionTemplate.h"

#include "CompartmentGrid.h"
#include "ChemCallbacks.h"

#include "Cylinder.h"
#include "Bead.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

void PolyPlusEndTemplate::addReaction(CCylinder* cc) {
    
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
            reactantSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(r)));
        
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
        for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
            productSpecies.push_back(m1->speciesBound(j));
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        
        if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+3>(species, _rate);
        else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+2>(species, _rate);
        else
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+1>(species, _rate);
        
        //callback
#ifdef REACTION_SIGNALING
        FilamentPolymerizationPlusEndCallback polyCallback(cc->getCylinder());
        ConnectionBlock rcb(rxn->connect(polyCallback,false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::POLYMERIZATIONPLUSEND);

        rxn->setRevNumber(_revnum);
        
    }
    
    //add extension callback reaction
    CMonomer* m = cc->getCMonomer(maxlength - 1);
    vector<Species*> reactantSpecies;
    vector<Species*> productSpecies;
    
    //loop through reactants, products. find all species
    auto r = _reactants[0];
    //FIRST REACTANT MUST BE BULK OR DIFFUSING
    if (getType(r) == SpeciesType::BULK)
        reactantSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(r)));
    
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
    for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
        productSpecies.push_back(m->speciesBound(j));
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    
    ReactionBase* rxn;
    if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+2>(species, _rate);
    else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+1>(species, _rate);
    else
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS>(species, _rate);
    
    //callbacks
#ifdef REACTION_SIGNALING
    short plusEndProduct = getInt(_products[1]);
    FilamentExtensionPlusEndCallback extCallback(cc->getCylinder(), plusEndProduct);
    ConnectionBlock rcb(rxn->connect(extCallback,false));
#endif
    
    cc->addInternalReaction(rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATIONPLUSEND);
    rxn->setRevNumber(_revnum);

}

void PolyMinusEndTemplate::addReaction(CCylinder* cc) {
    
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
            reactantSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(r)));
        
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
        for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
            productSpecies.push_back(m1->speciesBound(j));
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase *rxn;
        if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+3>(species, _rate);
        else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+2>(species, _rate);
        else
            rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+1>(species, _rate);
        
#ifdef REACTION_SIGNALING
        FilamentPolymerizationMinusEndCallback polyCallback(cc->getCylinder());
        ConnectionBlock rcb(rxn->connect(polyCallback,false));
#endif
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
        reactantSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(r)));
    
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
    for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
        productSpecies.push_back(m->speciesBound(j));
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    
    ReactionBase* rxn;
    if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+2>(species, _rate);
    else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS+1>(species, _rate);
    else
        rxn = new Reaction<POLYREACTANTS,POLYPRODUCTS>(species, _rate);
    
#ifdef REACTION_SIGNALING
    auto minusEndType = get<0>(_products[1]);
    FilamentExtensionMinusEndCallback extCallback(cc->getCylinder(), minusEndType);
    ConnectionBlock rcb(rxn->connect(extCallback,false));
#endif
    
    cc->addInternalReaction(rxn);
    rxn->setReactionType(ReactionType::POLYMERIZATIONMINUSEND);

    rxn->setRevNumber(_revnum);
}

void DepolyPlusEndTemplate::addReaction(CCylinder* cc) {
    
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
        for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
            reactantSpecies.push_back(m2->speciesBound(j));
        
        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        if( getType(p) == SpeciesType::BULK)
            productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
        
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS END
        p = _products[1];
        productSpecies.push_back(m2->speciesPlusEnd(getInt(p)));
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        
        if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
            rxn = new Reaction<DEPOLYREACTANTS+3,DEPOLYPRODUCTS>(species, _rate);
        else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
            rxn = new Reaction<DEPOLYREACTANTS+2,DEPOLYPRODUCTS>(species, _rate);
        else
            rxn = new Reaction<DEPOLYREACTANTS+1,DEPOLYPRODUCTS>(species, _rate);
        
#ifdef REACTION_SIGNALING
        FilamentDepolymerizationPlusEndCallback depolyCallback(cc->getCylinder());
        ConnectionBlock rcb(rxn->connect(depolyCallback,false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATIONPLUSEND);
        rxn->setRevNumber(_revnum);
    }
}

void DepolyMinusEndTemplate::addReaction(CCylinder* cc) {
    
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
        for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
            reactantSpecies.push_back(m2->speciesBound(j));
        
        //FIRST PRODUCT MUST BE BULK OR DIFFUSING
        auto p = _products[0];
        if(getType(p) == SpeciesType::BULK)
            productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
        
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        //SECOND PRODUCT SPECIES MUST BE MINUSEND
        p = _products[1];
        productSpecies.push_back(m2->speciesMinusEnd(getInt(p)));
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        
        if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
            rxn = new Reaction<DEPOLYREACTANTS+3,DEPOLYPRODUCTS>(species, _rate);
        else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
            rxn = new Reaction<DEPOLYREACTANTS+2,DEPOLYPRODUCTS>(species, _rate);
        else
            rxn = new Reaction<DEPOLYREACTANTS+1,DEPOLYPRODUCTS>(species, _rate);
        
#ifdef REACTION_SIGNALING
        FilamentDepolymerizationMinusEndCallback depolyCallback(cc->getCylinder());
        ConnectionBlock rcb(rxn->connect(depolyCallback,false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::DEPOLYMERIZATIONMINUSEND);

        rxn->setRevNumber(_revnum);
    }
}

void DepolyPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
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
    for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
        reactantSpecies.push_back(m2->speciesBound(j));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
    
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc2->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    //SECOND PRODUCT SPECIES MUST BE PLUS END
    p = _products[1];
    productSpecies.push_back(m2->speciesPlusEnd(getInt(p)));
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    
    ReactionBase* rxn;
    if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
        rxn = new Reaction<DEPOLYREACTANTS+3,DEPOLYPRODUCTS>(species, _rate);
    else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
        rxn = new Reaction<DEPOLYREACTANTS+2,DEPOLYPRODUCTS>(species, _rate);
    else
        rxn = new Reaction<DEPOLYREACTANTS+1,DEPOLYPRODUCTS>(species, _rate);
    
#ifdef REACTION_SIGNALING
    FilamentRetractionPlusEndCallback retCallback(cc1->getCylinder());
    ConnectionBlock rcb(rxn->connect(retCallback,false));
#endif
    
    cc2->addCrossCylinderReaction(cc1, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATIONPLUSEND);
    rxn->setRevNumber(_revnum);
}

void DepolyMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
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
    for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
        reactantSpecies.push_back(m2->speciesBound(j));
    
    //FIRST PRODUCT MUST BE BULK OR DIFFUSING
    auto p = _products[0];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
    
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    //SECOND PRODUCT SPECIES MUST BE MINUSEND
    p = _products[1];
    productSpecies.push_back(m2->speciesMinusEnd(getInt(p)));
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    
    ReactionBase* rxn;
    
    if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
        rxn = new Reaction<DEPOLYREACTANTS+3,DEPOLYPRODUCTS>(species, _rate);
    else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
        rxn = new Reaction<DEPOLYREACTANTS+2,DEPOLYPRODUCTS>(species, _rate);
    else
        rxn = new Reaction<DEPOLYREACTANTS+1,DEPOLYPRODUCTS>(species, _rate);
    
#ifdef REACTION_SIGNALING
    FilamentRetractionMinusEndCallback retCallback(cc1->getCylinder());
    ConnectionBlock rcb(rxn->connect(retCallback,false));
#endif
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::DEPOLYMERIZATIONMINUSEND);
    rxn->setRevNumber(_revnum);
}

void MotorWalkPTemplate::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
             it != SysParams::Chemistry().bindingSites[_filamentType].end() - 1; it++) {
        
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
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
        new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
        
        //callbacks
#ifdef REACTION_SIGNALING
        MotorWalkingCallback
        motorMoveCallback(cc->getCylinder(), site1, site2,
                          motorType, boundType, _ps);
        ConnectionBlock rcb(rxn->connect(motorMoveCallback, false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::MOTORWALKINGFORWARD);
    }
}

void MotorWalkPTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    CMonomer* m1 = cc1->getCMonomer(SysParams::Chemistry().bindingSites[_filamentType].back());
    CMonomer* m2 = cc2->getCMonomer(SysParams::Chemistry().bindingSites[_filamentType].front());
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
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn =
    new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
    
    //callbacks
#ifdef REACTION_SIGNALING
    MotorMovingCylinderCallback
    motorChangeCallback(cc1->getCylinder(), cc2->getCylinder(),
                        SysParams::Chemistry().bindingSites[_filamentType].back(),
                        SysParams::Chemistry().bindingSites[_filamentType].front(),
                        motorType, boundType, _ps);
    ConnectionBlock rcb(rxn->connect(motorChangeCallback, false));
#endif
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::MOTORWALKINGFORWARD);

}

void MotorWalkMTemplate::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = SysParams::Chemistry().bindingSites[_filamentType].end() - 1;
             it != SysParams::Chemistry().bindingSites[_filamentType].begin(); it--) {
        
        int site1 = *(it);
        int site2 = *(it-1);
        
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
        
        //Add the reaction. If it needs a callback then attach
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
        new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
        
        //callbacks
#ifdef REACTION_SIGNALING
        MotorWalkingCallback
        motorMoveCallback(cc->getCylinder(), site1, site2,
                          motorType, boundType, _ps);
        ConnectionBlock rcb(rxn->connect(motorMoveCallback, false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::MOTORWALKINGBACKWARD);
    }
}

void MotorWalkMTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
    CMonomer* m1 = cc2->getCMonomer(SysParams::Chemistry().bindingSites[_filamentType].front());
    CMonomer* m2 = cc1->getCMonomer(SysParams::Chemistry().bindingSites[_filamentType].back());
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
    
    //Add the reaction. If it needs a callback then attach
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn =
    new Reaction<MWALKINGREACTANTS, MWALKINGPRODUCTS>(species, _rate);
    
    //callbacks
#ifdef REACTION_SIGNALING
    MotorMovingCylinderCallback
    motorChangeCallback(cc2->getCylinder(), cc1->getCylinder(),
                        SysParams::Chemistry().bindingSites[_filamentType].front(),
                        SysParams::Chemistry().bindingSites[_filamentType].back(),
                        motorType, boundType, _ps);
    ConnectionBlock rcb(rxn->connect(motorChangeCallback, false));
#endif
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::MOTORWALKINGBACKWARD);
}

void AgingTemplate::addReaction(CCylinder* cc) {
    
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
        
        if(type == SpeciesType::FILAMENT) reactantSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)  reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND) reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //FIRST PRODUCT MUST BE FILAMENT, PLUS OR MINUS SPECIES
        auto p = _products[0];
        type = getType(p);
        speciesInt = getInt(p);
        
        if(type == SpeciesType::FILAMENT) productSpecies.push_back(m1->speciesFilament(speciesInt));
        else if(type == SpeciesType::PLUSEND)  productSpecies.push_back(m1->speciesPlusEnd(speciesInt));
        else if(type == SpeciesType::MINUSEND) productSpecies.push_back(m1->speciesMinusEnd(speciesInt));
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
        new Reaction<AGINGREACTANTS,AGINGPRODUCTS>(species, _rate);
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::AGING);
        rxn->setRevNumber(_revnum);
        
    }
}


void DestructionTemplate::addReaction(CCylinder* cc) {
    
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
            productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
        
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        p = _products[1];
        
        if(getType(p) == SpeciesType::BULK)
            productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
        
        else if(getType(p) == SpeciesType::DIFFUSING) {
            Compartment* c = cc->getCompartment();
            productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
        }
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* rxn =
        new Reaction<DESTRUCTIONREACTANTS,DESTRUCTIONPRODUCTS>(species, _rate);
        
#ifdef REACTION_SIGNALING
        FilamentDestructionCallback dcallback(cc->getCylinder(), _ps);
        ConnectionBlock rcb(rxn->connect(dcallback,false));
#endif
        
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::FILAMENTDESTRUCTION);
    }
}

void DestructionTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
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
        productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
    
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    p = _products[1];
    if(getType(p) == SpeciesType::BULK)
        productSpecies.push_back(_ps->getCompartmentGrid()->findSpeciesBulkByMolecule(getInt(p)));
    
    else if(getType(p) == SpeciesType::DIFFUSING) {
        Compartment* c = cc1->getCompartment();
        productSpecies.push_back(c->findSpeciesByMolecule(getInt(p)));
    }
    
    //Add the reaction
    vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* rxn =
    new Reaction<DESTRUCTIONREACTANTS,DESTRUCTIONPRODUCTS>(species, _rate);
    
#ifdef REACTION_SIGNALING
    FilamentDestructionCallback dcallback(cc1->getCylinder(), _ps);
    ConnectionBlock rcb(rxn->connect(dcallback,false));
#endif
    
    cc1->addCrossCylinderReaction(cc2, rxn);
    rxn->setReactionType(ReactionType::FILAMENTDESTRUCTION);
}

void SeveringTemplate::addReaction(CCylinder* cc) {
    
    //loop through all monomers
    for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
             it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
        
        int site = *(it);
        CMonomer* m = cc->getCMonomer(site);
        vector<Species*> reactantSpecies;
        
        //REACTANT MUST BE FILAMENT
        auto r = _reactants[0];
        reactantSpecies.push_back(m->speciesFilament(getInt(r)));
        
        //IMPLICITLY NEEDS AN EMPTY BOUND
        for(auto j : SysParams::Chemistry().bindingIndices[_filamentType])
            reactantSpecies.push_back(m->speciesBound(j));
        
        //Add the reaction
        vector<Species*> species = reactantSpecies;
        
        ReactionBase* rxn;
        if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 3)
            rxn = new Reaction<SEVERINGREACTANTS + 3,SEVERINGPRODUCTS>(species, _rate);
        else if(SysParams::Chemistry().bindingIndices[_filamentType].size() == 2)
            rxn = new Reaction<SEVERINGREACTANTS + 2,SEVERINGPRODUCTS>(species, _rate);
        else
            rxn = new Reaction<SEVERINGREACTANTS + 1,SEVERINGPRODUCTS>(species, _rate);
        
#ifdef REACTION_SIGNALING
        FilamentSeveringCallback scallback(cc->getCylinder());
        ConnectionBlock rcb(rxn->connect(scallback,false));
#endif
        cc->addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::SEVERING);
    }
}

SubSystem* FilamentReactionTemplate::_ps = 0;
