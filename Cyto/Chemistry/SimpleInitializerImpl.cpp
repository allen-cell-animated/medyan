//
//  SimpleInitializerImpl.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SimpleInitializerImpl.h"
#include "SystemParameters.h"
#include "ChemCallbacks.h"

void SimpleInitializerImpl::createFilamentReactionTemplates(ChemistrySpeciesAndReactions& chemSR) {
    
    ///set up reaction templates
    for(auto &r: chemSR.filamentReactions) {
        
        std::vector<std::tuple<int, SpeciesType>> reactantTemplate;
        std::vector<std::tuple<int, SpeciesType>> productTemplate;
        FilamentReactionDirection d;
        ReactionType type;
        
        for(std::string& reactant : std::get<0>(r)) {
            
            ///read strings, and look up type
            if(reactant.find("BULK") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find_if(chemSR.speciesBulk.begin(), chemSR.speciesBulk.end(),
                                       [name](std::tuple<std::string, int> element) { return std::get<0>(element) == name ? true : false; });
                                           
                if(it == chemSR.speciesBulk.end()) {
                    std::cout << "A bulk species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(reactant.find("DIFFUSING") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find_if(chemSR.speciesDiffusing.begin(), chemSR.speciesDiffusing.end(),
                                       [name](std::tuple<std::string, int, double> element) { return std::get<0>(element) == name ? true : false; });
                if(it == chemSR.speciesDiffusing.end()) {
                    std::cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
            }
            
            else if(reactant.find("FILAMENT") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesFilament.begin(), _speciesFilament.end(), name);
                int position = 0;
                
                if(it != _speciesFilament.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesFilament.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    std::cout << "A filament species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(reactant.find("BOUND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesBound.begin(), _speciesBound.end(), name);
                int position = 0;
                
                if(it != _speciesBound.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesBound.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::BOUND));
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(reactant.find("LINKER") != std::string::npos) {
                
                type = ReactionType::LINKERUNBINDING;
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesLinker.begin(), _speciesLinker.end(), name);
                int position = 0;
                
                if(it != _speciesLinker.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesLinker.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::LINKER));
                }
                else {
                    std::cout << "A linker species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(reactant.find("MOTOR") != std::string::npos) {
                
                type = ReactionType::MOTORUNBINDING;
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesMotor.begin(), _speciesMotor.end(), name);
                int position = 0;
                
                if(it != _speciesMotor.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesMotor.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
                }
                else {
                    std::cout << "A motor species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(reactant.find("PLUSEND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
                int position = 0;
                
                if(it != _speciesPlusEnd.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesPlusEnd.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                    
                    ///see what position this is (N or N+1)
                    if(reactant.find("N+1") != std::string::npos) d = FilamentReactionDirection::BACKWARD;
                    else d = FilamentReactionDirection::FORWARD;
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(reactant.find("MINUSEND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
                int position = 0;
                
                if(it != _speciesMinusEnd.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesMinusEnd.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                    
                    ///see what position this is (N or N+1)
                    if(reactant.find("N+1") != std::string::npos) d = FilamentReactionDirection::BACKWARD;
                    else d = FilamentReactionDirection::FORWARD;
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        for(std::string& product : std::get<1>(r)) {
            
            ///read strings, and look up type
            if(product.find("BULK") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find_if(chemSR.speciesBulk.begin(), chemSR.speciesBulk.end(),
                                       [name](std::tuple<std::string, int> element) { return std::get<0>(element) == name ? true : false; });
                if(it == chemSR.speciesBulk.end()) {
                    std::cout << "A bulk species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(product.find("DIFFUSING") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find_if(chemSR.speciesDiffusing.begin(), chemSR.speciesDiffusing.end(),
                                       [name](std::tuple<std::string, int, double> element) { return std::get<0>(element) == name ? true : false; });
                if(it == chemSR.speciesDiffusing.end()) {
                    std::cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
            }
            
            else if(product.find("FILAMENT") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesFilament.begin(), _speciesFilament.end(), name);
                int position = 0;
                
                if(it != _speciesFilament.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesFilament.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    std::cout << "A filament species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(product.find("BOUND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesBound.begin(), _speciesBound.end(), name);
                int position = 0;
                
                if(it != _speciesBound.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesBound.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::BOUND));
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(product.find("PLUSEND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
                int position = 0;
                
                if(it != _speciesPlusEnd.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesPlusEnd.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                    
                    ///see what position this is (N or N+1)
                    if(product.find("N+1") != std::string::npos && d != FilamentReactionDirection::FORWARD) {
                        std::cout << "A filament reaction involving polymerization/depolymerization is invalid. Exiting." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    if(product.find("N+1") != std::string::npos) { type = ReactionType::POLYMERIZATION; }
                    else { type = type = ReactionType::DEPOLYMERIZATION; }
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            else if(product.find("MINUSEND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
                int position = 0;
                
                if(it != _speciesMinusEnd.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesMinusEnd.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                    
                    ///see what position this is (N or N+1)
                    if(product.find("N+1") != std::string::npos && d != FilamentReactionDirection::FORWARD) {
                        std::cout << "A filament reaction involving polymerization/depolymerization is invalid. Exiting." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    if(product.find("N+1") != std::string::npos) { type = ReactionType::DEPOLYMERIZATION; }
                    else { type = ReactionType::POLYMERIZATION; }
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        ///add reaction template
        if(type == ReactionType::POLYMERIZATION){
            
            if(d == FilamentReactionDirection::FORWARD)
                _reactionFilamentTemplates.emplace_back(
                    new PolymerizationPlusEndTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
            else
                _reactionFilamentTemplates.emplace_back(
                    new PolymerizationMinusEndTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
        }
        
        else if(type == ReactionType::DEPOLYMERIZATION){
           
            if(d == FilamentReactionDirection::BACKWARD)
                _reactionFilamentTemplates.emplace_back(
                    new DepolymerizationPlusEndTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
            else
                _reactionFilamentTemplates.emplace_back(
                    new DepolymerizationMinusEndTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
        }
        else if(type == ReactionType::LINKERUNBINDING) {
            _reactionFilamentTemplates.emplace_back(new LinkerUnbindingTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
        }
        else if(type == ReactionType::MOTORUNBINDING) {
            _reactionFilamentTemplates.emplace_back(new MotorUnbindingTemplate(reactantTemplate, productTemplate, std::get<2>(r)));
        }
    }
}

void SimpleInitializerImpl::createCrossFilamentReactionTemplates(ChemistrySpeciesAndReactions& chemSR) {
    
    std::vector<std::tuple<int, SpeciesType>> reactantTemplate;
    std::vector<std::tuple<int, SpeciesType>> productTemplate;
    ReactionType type;
    
    for(auto &r: chemSR.crossFilamentReactions) {
    
        for(std::string& reactant : std::get<0>(r)) {
            
            ///read strings, and look up type
            if(reactant.find("BULK") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find_if(chemSR.speciesBulk.begin(), chemSR.speciesBulk.end(),
                                       [name](std::tuple<std::string, int> element) { return std::get<0>(element) == name ? true : false; });
                
                if(it == chemSR.speciesBulk.end()) {
                    std::cout << "A bulk species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(reactant.find("DIFFUSING") != std::string::npos) {
                
                ///Look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find_if(chemSR.speciesDiffusing.begin(), chemSR.speciesDiffusing.end(),
                                       [name](std::tuple<std::string, int, double> element) { return std::get<0>(element) == name ? true : false; });
                if(it == chemSR.speciesDiffusing.end()) {
                    std::cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(std::tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
            }
            
            else if(reactant.find("BOUND") != std::string::npos) {
                
                ///look up species, make sure in list
                std::string name = reactant.substr(0, reactant.find(":"));
                auto it = std::find(_speciesBound.begin(), _speciesBound.end(), name);
                int position = 0;
                
                if(it != _speciesBound.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesBound.begin(), it);
                    reactantTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::BOUND));
                }
                else {
                    std::cout << "A bound species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        for(std::string& product : std::get<1>(r)) {
            
            if(product.find("LINKER") != std::string::npos) {
                
                type = ReactionType::LINKERBINDING;
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesLinker.begin(), _speciesLinker.end(), name);
                int position = 0;
                
                if(it != _speciesLinker.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesLinker.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::LINKER));
                }
                else {
                    std::cout << "A linker species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            if(product.find("MOTOR") != std::string::npos) {
                
                type = ReactionType::MOTORBINDING;
                
                ///look up species, make sure in list
                std::string name = product.substr(0, product.find(":"));
                auto it = std::find(_speciesMotor.begin(), _speciesMotor.end(), name);
                int position = 0;
                
                if(it != _speciesMotor.end()) {
                    
                    ///get position of iterator
                    position = std::distance(_speciesMotor.begin(), it);
                    productTemplate.push_back(std::tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
                }
                else {
                    std::cout << "A motor species that was included in a reaction was not initialized. Exiting." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    
        if(type == ReactionType::LINKERBINDING)
            _reactionFilamentTemplates.emplace_back(
                new LinkerBindingTemplate(reactantTemplate, productTemplate, std::get<2>(r), std::get<3>(r), std::get<4>(r)));
        
        else if(type == ReactionType::MOTORBINDING)
            _reactionFilamentTemplates.emplace_back(
                new MotorBindingTemplate(reactantTemplate, productTemplate, std::get<2>(r), std::get<3>(r), std::get<4>(r)));
    }
}


void SimpleInitializerImpl::initialize(ChemistrySpeciesAndReactions& chemSR) {
    
    ///Copy all species from chemSR struct
    _speciesFilament = chemSR.speciesFilament; _speciesPlusEnd = chemSR.speciesPlusEnd; _speciesMinusEnd = chemSR.speciesMinusEnd;
    
    ///add bound, linkers and motors
    _speciesBound = chemSR.speciesBound; _speciesLinker = chemSR.speciesLinker; _speciesMotor = chemSR.speciesMotor;
    
    ///Setup all species diffusing and bulk
    Compartment& cProto = CompartmentGrid::Instance(compartmentGridKey())->getProtoCompartment();
    
    for(auto &sb : chemSR.speciesBulk)
        CompartmentGrid::Instance(compartmentGridKey())->addSpeciesBulk(std::get<0>(sb), std::get<1>(sb));
    for(auto &sd : chemSR.speciesDiffusing)
        cProto.addSpeciesDiffusing(std::get<0>(sd), std::get<1>(sd), std::get<2>(sd));
    
    ///initialize all compartments with species diffusing
    for(auto &c : CompartmentGrid::Instance(compartmentGridKey())->children()) {
        Compartment *C = static_cast<Compartment*>(c.get());
        *C = cProto;
    }
    ///activate all compartments for diffusion, set up diffusion reactions
    CompartmentGrid::Instance(compartmentGridKey())->activateAll();
    for(auto &c : CompartmentGrid::Instance(compartmentGridKey())->children()) {
        Compartment *C = static_cast<Compartment*>(c.get());
        C->generateAllDiffusionReactions();
    }
    
    ///add reactions to chemsim
    CompartmentGrid::Instance(compartmentGridKey())->addChemSimReactions();
    
    ///create filament reaction templates
    createFilamentReactionTemplates(chemSR);
    ///create cross filament reaction templates
    createCrossFilamentReactionTemplates(chemSR);
}

CCylinder* SimpleInitializerImpl::createCCylinder(Filament *pf, Compartment* c,
                                                  bool extensionFront, bool extensionBack)
{
    CCylinder* cc = new CCylinder(c);
    
    ///maxlength is same length as mcylinder
    int maxlength = cc->size();
    
    ///add monomers to cylinder
    for(int i = 0; i < maxlength; i++) {
        
        CMonomer* m = new CMonomer();
        for(auto &f : _speciesFilament) {
            SpeciesFilament* sf =
                c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName(f), 0, 1);
            m->addSpeciesFilament(sf);
        }
        for (auto &p : _speciesPlusEnd) {
            SpeciesPlusEnd* sp =
                c->addSpeciesPlusEnd(SpeciesNamesDB::Instance()->generateUniqueName(p), 0, 1);
            m->addSpeciesPlusEnd(sp);
        }
        for (auto &mi : _speciesMinusEnd) {
            SpeciesMinusEnd* smi =
                c->addSpeciesMinusEnd(SpeciesNamesDB::Instance()->generateUniqueName(mi), 0, 1);
            m->addSpeciesMinusEnd(smi);
        }
        
        for (auto &b : _speciesBound) {
            SpeciesBound* sb =
                c->addSpeciesBound(SpeciesNamesDB::Instance()->generateUniqueName(b), 0, 1);
            m->addSpeciesBound(sb);
        }
        for (auto &l : _speciesLinker) {
            SpeciesLinker* sl =
                c->addSpeciesLinker(SpeciesNamesDB::Instance()->generateUniqueName(l), 0, 1);
            m->addSpeciesLinker(sl);
        }
        for (auto &mo : _speciesMotor) {
            SpeciesMotor* sm =
                c->addSpeciesMotor(SpeciesNamesDB::Instance()->generateUniqueName(mo), 0, 1);
            m->addSpeciesMotor(sm);
        }
        
        cc->addCMonomer(m);
    }
    
    ///Add all reaction templates to this cylinder
    for(auto &r : _reactionFilamentTemplates) { r->addReaction(cc, pf); }
    
    ///get last ccylinder
    CCylinder* lastcc = nullptr;
 
    ///extension of front
    if(extensionFront) {
        lastcc = pf->getCylinderVector().back()->getCCylinder();
        for(auto &r : _reactionFilamentTemplates) r->addReaction(lastcc, cc, pf);
    }
    ///extension of back
    else if(extensionBack) {
        lastcc = pf->getCylinderVector().front()->getCCylinder();
        for(auto &r : _reactionFilamentTemplates) r->addReaction(cc, lastcc, pf);
    }

    ///Base case, initialization
    else {
        ///Check if this is the first cylinder
        if(pf->getCylinderVector().size() != 0) {
            
            ///remove plus end from last, add to this.
            lastcc = pf->getCylinderVector().back()->getCCylinder();
            CMonomer* m1 = lastcc->getCMonomer(lastcc->size() - 2);
            m1->speciesPlusEnd(0)->getRSpecies().setN(0);
            
            CMonomer* m2 = cc->getCMonomer(cc->size() - 2);
            m2->speciesPlusEnd(0)->getRSpecies().setN(1);
            m2->speciesBound(0)->getRSpecies().setN(1);
            
            ///fill last cylinder with default filament value
            for(int i = lastcc->size() - 2; i < lastcc->size(); i++) {
                lastcc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                lastcc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
                
            }
            ///fill new cylinder with default filament value
            for(int i = 0; i < cc->size() - 2; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
            
            for(auto &r : _reactionFilamentTemplates) r->addReaction(lastcc, cc, pf);
            
        }
        ///this is first one
        else {
            //set back and front
            CMonomer* m1 = cc->getCMonomer(cc->size() - 2);
            m1->speciesPlusEnd(0)->getRSpecies().setN(1);
            m1->speciesBound(0)->getRSpecies().setN(1);
            
            CMonomer* m2 = cc->getCMonomer(1);
            m2->speciesMinusEnd(0)->getRSpecies().setN(1);
            m2->speciesBound(0)->getRSpecies().setN(1);
            
            ///fill with default filament value
            for(int i = 2; i < cc->size() - 2; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
        }
    }
    
    //update all reactions added
    cc->updateReactions();
    
    ///cc->printCCylinder();
    //std::cout <<std::endl;

    ///clean up and return
    return cc;
}

