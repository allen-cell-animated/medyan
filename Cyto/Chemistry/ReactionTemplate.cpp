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
        for(auto &r : _reactants) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                
                case SpeciesType::BULK: {
                    reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                                         findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::PLUSEND: {
                    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
                    break;
                }
                default: {}
            }
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::FILAMENT: {
                    productSpecies.push_back(m1->speciesFilament(speciesInt));
                    break;
                }
                    
                case SpeciesType::BOUND: {
                    productSpecies.push_back(m2->speciesBound(speciesInt));
                    break;
                }
                    
                case SpeciesType::PLUSEND: {
                    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));
                    break;
                }
                    
                default: {}
            }
        }
        
        FilamentExtensionFrontCallback extCallback(pf);
        FilamentPolymerizationFrontCallback polyCallback(pf);
        
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<2, 3>(species, _rate);

        if(i == maxlength - 2)
            boost::signals2::shared_connection_block rcb(r->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block rcb(r->connect(polyCallback,false));
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::POLYMERIZATION);
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
        for(auto &r : _reactants) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::BULK: {
                    reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                              findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::MINUSEND: {
                    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
                    break;
                }
                default: {}
            }
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::FILAMENT: {
                    productSpecies.push_back(m1->speciesFilament(speciesInt));
                    break;
                }
                    
                case SpeciesType::BOUND: {
                    productSpecies.push_back(m2->speciesBound(speciesInt));
                    break;
                }
                    
                case SpeciesType::MINUSEND: {
                    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
                    break;
                }
                    
                default: {}
            }
        }
        
        FilamentExtensionBackCallback extCallback(pf);
        FilamentPolymerizationBackCallback polyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<2, 3>(species, _rate);
        
        if(i == 1)
            boost::signals2::shared_connection_block rcb(r->connect(extCallback,false));
        else
            boost::signals2::shared_connection_block rcb(r->connect(polyCallback,false));
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::POLYMERIZATION);
    }

}



void PolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    for(auto &r : _reactants) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::BULK: {
                reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                          findSpeciesBulkByMolecule(speciesInt));
                break;
            }
            case SpeciesType::DIFFUSING: {
                Compartment* c = cc1->getCompartment();
                reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                break;
            }
            case SpeciesType::PLUSEND: {
                reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
                break;
            }
            default: {}
        }
    }
    
    for(auto &r: _products) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::FILAMENT: {
                productSpecies.push_back(m1->speciesFilament(speciesInt));
                break;
            }
                
            case SpeciesType::BOUND: {
                productSpecies.push_back(m2->speciesBound(speciesInt));
                break;
            }
                
            case SpeciesType::PLUSEND: {
                productSpecies.push_back(m2->speciesPlusEnd(speciesInt));
                break;
            }
                
            default: {}
        }
    }
    
    FilamentPolymerizationFrontCallback polyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(r->connect(polyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, r);
    r->setReactionType(ReactionType::POLYMERIZATION);
}

void PolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {
    
    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    for(auto &r : _reactants) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::BULK: {
                reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                          findSpeciesBulkByMolecule(speciesInt));
                break;
            }
            case SpeciesType::DIFFUSING: {
                Compartment* c = cc2->getCompartment();
                reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                break;
            }
            case SpeciesType::MINUSEND: {
                reactantSpecies.push_back(m2->speciesMinusEnd(speciesInt));
                break;
            }
            default: {}
        }
    }
    
    for(auto &r: _products) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::FILAMENT: {
                productSpecies.push_back(m2->speciesFilament(speciesInt));
                break;
            }
                
            case SpeciesType::BOUND: {
                productSpecies.push_back(m1->speciesBound(speciesInt));
                break;
            }
                
            case SpeciesType::MINUSEND: {
                productSpecies.push_back(m1->speciesMinusEnd(speciesInt));
                break;
            }
                
            default: {}
        }
    }
    
    FilamentPolymerizationBackCallback polyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<2, 3>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(r->connect(polyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, r);
    r->setReactionType(ReactionType::POLYMERIZATION);
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
        for(auto &r : _reactants) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::FILAMENT: {
                    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
                    break;
                }
                    
                case SpeciesType::BOUND: {
                    reactantSpecies.push_back(m1->speciesBound(speciesInt));
                    break;
                }
                    
                case SpeciesType::PLUSEND: {
                    reactantSpecies.push_back(m1->speciesPlusEnd(speciesInt));
                    break;
                }
                default: {}
            }
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::BULK: {
                    productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                              findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::PLUSEND: {
                    productSpecies.push_back(m2->speciesPlusEnd(speciesInt));
                    break;
                }
                default: {}
            }
        }
        
        FilamentRetractionFrontCallback retCallback(pf);
        FilamentDepolymerizationFrontCallback depolyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<3, 2>(species, _rate);
        
        if(i == maxlength - 1)
            boost::signals2::shared_connection_block rcb(r->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block rcb(r->connect(depolyCallback,false));
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::DEPOLYMERIZATION);
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
        for(auto &r : _reactants) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::FILAMENT: {
                    reactantSpecies.push_back(m2->speciesFilament(speciesInt));
                    break;
                }
                    
                case SpeciesType::BOUND: {
                    reactantSpecies.push_back(m1->speciesBound(speciesInt));
                    break;
                }
                    
                case SpeciesType::MINUSEND: {
                    reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
                    break;
                }
                default: {}
            }
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::BULK: {
                    productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                             findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::MINUSEND: {
                    productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
                    break;
                }
                    
                default: {}
            }
        }
        
        FilamentRetractionBackCallback retCallback(pf);
        FilamentDepolymerizationBackCallback depolyCallback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<3, 2>(species, _rate);
        
        if(i == 0)
            boost::signals2::shared_connection_block rcb(r->connect(retCallback,false));
        else
            boost::signals2::shared_connection_block rcb(r->connect(depolyCallback,false));
            
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::DEPOLYMERIZATION);
    }
}

void DepolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {
    

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    for(auto &r : _reactants) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::FILAMENT: {
                reactantSpecies.push_back(m1->speciesFilament(speciesInt));
                break;
            }
                
            case SpeciesType::BOUND: {
                reactantSpecies.push_back(m2->speciesBound(speciesInt));
                break;
            }
                
            case SpeciesType::PLUSEND: {
                reactantSpecies.push_back(m2->speciesPlusEnd(speciesInt));
                break;
            }
            default: {}
        }
    }
    
    for(auto &r: _products) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
                
            case SpeciesType::BULK: {
                productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                          findSpeciesBulkByMolecule(speciesInt));
                break;
            }
            case SpeciesType::DIFFUSING: {
                Compartment* c = cc2->getCompartment();
                productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                break;
            }
            case SpeciesType::PLUSEND: {
                productSpecies.push_back(m1->speciesPlusEnd(speciesInt));
                break;
            }
                
            default: {}
        }
    }
    FilamentDepolymerizationFrontCallback depolyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(r->connect(depolyCallback,false));
    
    cc2->addCrossCylinderReaction(cc1, r);
    r->setReactionType(ReactionType::DEPOLYMERIZATION);

}

void DepolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {

    CMonomer* m1 = cc1->getCMonomer(cc1->size() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    std::vector<Species*> reactantSpecies;
    std::vector<Species*> productSpecies;
    
    ///loop through reactants, products. find all species
    for(auto &r : _reactants) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
            case SpeciesType::FILAMENT: {
                reactantSpecies.push_back(m2->speciesFilament(speciesInt));
                break;
            }
                
            case SpeciesType::BOUND: {
                reactantSpecies.push_back(m1->speciesBound(speciesInt));
                break;
            }
                
            case SpeciesType::MINUSEND: {
                reactantSpecies.push_back(m1->speciesMinusEnd(speciesInt));
                break;
            }
            default: {}
        }
    }
    
    for(auto &r: _products) {
        
        SpeciesType type = std::get<1>(r);
        int speciesInt = std::get<0>(r);
        
        switch(type) {
                
                
            case SpeciesType::BULK: {
                productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                         findSpeciesBulkByMolecule(speciesInt));
                break;
            }
            case SpeciesType::DIFFUSING: {
                Compartment* c = cc1->getCompartment();
                productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                break;
            }
            case SpeciesType::MINUSEND: {
                productSpecies.push_back(m2->speciesMinusEnd(speciesInt));
                break;
            }
                
            default: {}
        }
    }
    
    FilamentDepolymerizationBackCallback depolyCallback(pf);
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<3, 2>(species, _rate);
    
    boost::signals2::shared_connection_block rcb(r->connect(depolyCallback,false));
    
    cc1->addCrossCylinderReaction(cc2, r);
    r->setReactionType(ReactionType::DEPOLYMERIZATION);
}

void LinkerUnbindingTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        for(auto &r : _reactants) {
            
            int speciesInt = std::get<0>(r);

            reactantSpecies.push_back(m1->speciesLinker(speciesInt));
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::BULK: {
                    productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                             findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::BOUND: {
                    productSpecies.push_back(m1->speciesBound(speciesInt));
                    break;
                }
                    
                default: {}
            }
        }
        LinkerUnbindingCallback unbindingCallback(cc);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<1, 2>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(r->connect(unbindingCallback,false));
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::LINKERUNBINDING);
    }   
}

void MotorUnbindingTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();
    
    ///loop through all monomers
    for(int i = 0; i < maxlength - 1; i++) {
        
        CMonomer* m1 = cc->getCMonomer(i);
        std::vector<Species*> reactantSpecies;
        std::vector<Species*> productSpecies;
        
        ///loop through reactants, products. find all species
        for(auto &r : _reactants) {
            
            int speciesInt = std::get<0>(r);
            
            reactantSpecies.push_back(m1->speciesMotor(speciesInt));
        }
        
        for(auto &r: _products) {
            
            SpeciesType type = std::get<1>(r);
            int speciesInt = std::get<0>(r);
            
            switch(type) {
                    
                case SpeciesType::BULK: {
                    productSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                             findSpeciesBulkByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::DIFFUSING: {
                    Compartment* c = cc->getCompartment();
                    productSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                    break;
                }
                case SpeciesType::BOUND: {
                    productSpecies.push_back(m1->speciesBound(speciesInt));
                    break;
                }
                    
                default: {}
            }
        }
        MotorUnbindingCallback unbindingCallback(cc);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<1, 2>(species, _rate);
        
        boost::signals2::shared_connection_block rcb(r->connect(unbindingCallback,false));
        
        cc->addInternalReaction(r);
        r->setReactionType(ReactionType::MOTORUNBINDING);
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
    
    if(_rMin <= dist && dist <= _rMax) {
        ///loop through all monomers of filament
        int maxlength1 = cc1->size();
        int maxlength2 = cc2->size();
        
        ///loop through all monomers
        for(int i = 0; i < maxlength1 - 1; i++) {
            
            CMonomer* m1 = cc1->getCMonomer(i);
            
            for(int j = 0; j < maxlength2 - 1; j++) {
            
                CMonomer* m2 = cc2->getCMonomer(j);
                std::vector<Species*> reactantSpecies;
                std::vector<Species*> productSpecies;
                
                int firstBound = 0;
                ///loop through reactants, products. find all species
                for(auto &r : _reactants) {
                    
                    SpeciesType type = std::get<1>(r);
                    int speciesInt = std::get<0>(r);
                    
                    switch(type) {
                            
                        case SpeciesType::BULK: {
                            reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                                      findSpeciesBulkByMolecule(speciesInt));
                            break;
                        }
                        case SpeciesType::DIFFUSING: {
                            Compartment* c = cc1->getCompartment();
                            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                            break;
                        }
                        case SpeciesType::BOUND: {
                            
                            if(firstBound == 0) reactantSpecies.push_back(m1->speciesBound(speciesInt));
                            else reactantSpecies.push_back(m2->speciesBound(speciesInt));
                            firstBound++;
                            break;
                        }
                        default: {}
                    }
                }
                firstBound = 0;
                for(auto &r: _products) {

                    int speciesInt = std::get<0>(r);

                    if(firstBound == 0) productSpecies.push_back(m1->speciesBound(speciesInt));
                    else productSpecies.push_back(m2->speciesBound(speciesInt));
                }
       
                ///set up callbacks
                LinkerBindingCallback lcallback(int linkerNumber, cc1, cc2, system);
                
                ///Add the reaction. If it needs a callback then attach
                std::vector<Species*> species = reactantSpecies;
                species.insert(species.end(), productSpecies.begin(), productSpecies.end());
                ReactionBase* r = new Reaction<3, 2>(species, _rate);
                
                boost::signals2::shared_connection_block rcb(r->connect(lcallback,false));
                
                cc1->addCrossCylinderReaction(cc2, r);
                r->setReactionType(ReactionType::LINKERBINDING);
            }
        }
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
    
    if(_rMin <= dist && dist <= _rMax) {
        ///loop through all monomers of filament
        int maxlength1 = cc1->size();
        int maxlength2 = cc2->size();
        
        ///loop through all monomers
        for(int i = 0; i < maxlength1 - 1; i++) {
            
            CMonomer* m1 = cc1->getCMonomer(i);
            
            for(int j = 0; j < maxlength2 - 1; j++) {
                
                CMonomer* m2 = cc2->getCMonomer(j);
                std::vector<Species*> reactantSpecies;
                std::vector<Species*> productSpecies;
                
                int firstBound = 0;
                ///loop through reactants, products. find all species
                for(auto &r : _reactants) {
                    
                    SpeciesType type = std::get<1>(r);
                    int speciesInt = std::get<0>(r);
                    
                    switch(type) {
                            
                        case SpeciesType::BULK: {
                            reactantSpecies.push_back(&CompartmentGrid::Instance(CompartmentGridKey())->
                                                      findSpeciesBulkByMolecule(speciesInt));
                            break;
                        }
                        case SpeciesType::DIFFUSING: {
                            Compartment* c = cc1->getCompartment();
                            reactantSpecies.push_back(c->findSpeciesByMolecule(speciesInt));
                            break;
                        }
                        case SpeciesType::BOUND: {
                            
                            if(firstBound == 0) reactantSpecies.push_back(m1->speciesBound(speciesInt));
                            else reactantSpecies.push_back(m2->speciesBound(speciesInt));
                            firstBound++;
                            break;
                        }
                        default: {}
                    }
                }
                firstBound = 0;
                for(auto &r: _products) {
                    
                    int speciesInt = std::get<0>(r);
                    
                    if(firstBound == 0) productSpecies.push_back(m1->speciesBound(speciesInt));
                    else productSpecies.push_back(m2->speciesBound(speciesInt));
                }
                
                ///set up callbacks
                MotorBindingCallback lcallback(int linkerNumber, cc1, cc2, system);
                
                ///Add the reaction. If it needs a callback then attach
                std::vector<Species*> species = reactantSpecies;
                species.insert(species.end(), productSpecies.begin(), productSpecies.end());
                ReactionBase* r = new Reaction<3, 2>(species, _rate);
                
                boost::signals2::shared_connection_block rcb(r->connect(lcallback,false));
                
                cc1->addCrossCylinderReaction(cc2, r);
                r->setReactionType(ReactionType::MOTORBINDING);
            }
        }
    }
}



