//
//  ReactionTemplate.cpp
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ReactionTemplate.h"
#include "SystemParameters.h"
#include "CompartmentContainer.h"
#include "ChemCallbacks.h"

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
        
        FilamentExtensionFrontCallback callback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<2, 3>(species, _rate);

        if(i == maxlength - 2)
            boost::signals2::shared_connection_block rcb1(r->connect(callback,false));
        
        cc->addReaction(r);
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
        
        FilamentExtensionBackCallback callback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<2, 3>(species, _rate);
        
        if(i == 1)
            boost::signals2::shared_connection_block rcb1(r->connect(callback,false));
        
        cc->addReaction(r);
    }

}



void PolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {

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
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<2, 3>(species, _rate);
    
    cc1->addFrontReaction(r, true);
    cc2->addBackReaction(r);

}

void PolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    
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
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<2, 3>(species, _rate);
    
    cc2->addBackReaction(r, true);
    cc1->addFrontReaction(r);


}


void DepolymerizationPlusEndTemplate::addReaction(CCylinder* cc, Filament* pf) {
    
    ///loop through all monomers of filament
    int maxlength = cc->size();

    ///loop through all monomers
    for(int i = maxlength - 1; i > 0; i--) {
        
        std::cout << i << std::endl;
        
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
        
        FilamentRetractionFrontCallback callback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<3, 2>(species, _rate);
        
        if(i == maxlength - 1)
            //boost::signals2::shared_connection_block rcb1(r->connect(callback,false));
        
        cc->addReaction(r);
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
        
        FilamentRetractionBackCallback callback(pf);
        
        ///Add the reaction. If it needs a callback then attach
        std::vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        ReactionBase* r = new Reaction<3, 2>(species, _rate);
        
        if(i == 0)
            boost::signals2::shared_connection_block rcb1(r->connect(callback,false));
        
        cc->addReaction(r);
    }

}



void DepolymerizationPlusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {
    

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
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<3, 2>(species, _rate);
    
    cc2->addBackReaction(r, true);
    cc1->addFrontReaction(r);


}

void DepolymerizationMinusEndTemplate::addReaction(CCylinder* cc1, CCylinder* cc2) {

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
    
    ///Add the reaction. If it needs a callback then attach
    std::vector<Species*> species = reactantSpecies;
    species.insert(species.end(), productSpecies.begin(), productSpecies.end());
    ReactionBase* r = new Reaction<3, 2>(species, _rate);
    
    cc1->addFrontReaction(r, true);
    cc2->addBackReaction(r);

}


