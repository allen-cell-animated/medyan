//
//  SimpleManagerImpl.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "SimpleManagerImpl.h"

#include "ChemCallbacks.h"
#include "ReactionManager.h"
#include "Parser.h"
#include "Bead.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

void SimpleManagerImpl::genIFRxnManagers(ChemistryData& chem) {
    
    ///set up reaction templates
    for(auto &r: chem.polymerizationReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        FilamentReactionDirection d;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        ///read strings, and look up type
        
        ///Checks on number of reactants, products
        if(reactants.size() != 2 || products.size() != 3) {
            cout << "Invalid polymerization reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST SPECIES MUST BE BULK OR DIFFUSING
        auto reactant = reactants[0];
        if(reactant.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
                                       
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "First species listed in a polymerization reaction must be either bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND SPECIES MUST BE PLUS OR MINUS END
        reactant = reactants[1];
        if(reactant.find("PLUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                
                d = FilamentReactionDirection::FORWARD;
            }
            else {
                cout << "A plus end species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(reactant.find("MINUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                
                d = FilamentReactionDirection::BACKWARD;
            }
            else {
                cout << "A minus species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            
            cout << "Second species listed in a polymerization reaction must be either plusend or minusend. Exiting." << endl;
            exit(EXIT_FAILURE);
            
        }
        
        ///FIRST PRODUCT SPECIES MUST BE FILAMENT SPECIES
        auto product = products[0];
        if(product.find("FILAMENT") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                ///get position of iterator
                position = distance(_speciesFilament.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
            }
            else {
                cout << "A filament species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Third species listed in a polymerization reaction must be filament. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND PRODUCT SPECIES MUST BE BOUND
        product = products[1];
        ///read strings, and look up type
        if(product.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Fourth species listed in a polymerization reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
            
        }
        
        ///THIRD PRODUCT SPECIES MUST BE PLUS OR MINUS END
        product = products[2];
        ///read strings, and look up type
        if(product.find("PLUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
            }
            else {
                cout << "A plusend species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(product.find("MINUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
            }
            else {
                cout << "A plusend species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Fifth species listed in a polymerization reaction must be either plusend or minusend. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///Add polymerization managers
        if(d == FilamentReactionDirection::FORWARD)
            _IFRxnManagers.emplace_back(new PolyPlusEndManager(reactantTemplate, productTemplate, get<2>(r)));
        else
            _IFRxnManagers.emplace_back(new PolyMinusEndManager(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    
    ///set up reaction templates
    for(auto &r: chem.depolymerizationReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        FilamentReactionDirection d;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        ///read strings, and look up type
        
        ///Checks on number of reactants, products
        if(reactants.size() != 3 || products.size() != 2) {
            cout << "Invalid depolymerization reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        ///FIRST REACTANT SPECIES MUST BE FILAMENT SPECIES
        auto reactant = reactants[0];
        if(reactant.find("FILAMENT") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                ///get position of iterator
                position = distance(_speciesFilament.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
            }
            else {
                cout << "A filament species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "First species listed in a depolymerization reaction must be filament. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND REACTANT SPECIES MUST BE BOUND
        reactant = reactants[1];
        ///read strings, and look up type
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Second species listed in a depolymerization reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
            
        }
        
        ///THIRD REACTANT SPECIES MUST BE PLUS OR MINUS END
        reactant = reactants[2];
        ///read strings, and look up type
        if(reactant.find("PLUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                
                d = FilamentReactionDirection::BACKWARD;
            }
            else {
                cout << "A plusend species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(reactant.find("MINUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                d = FilamentReactionDirection::FORWARD;
            }
            else {
                cout << "A minusend species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Third species listed in a depolymerization reaction must be either plusend or minusend. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        
        ///FIRST PRODUCT SPECIES MUST BE BULK OR DIFFUSING
        auto product = products[0];
        if(product.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(product.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Fourth species listed in a depolymerization reaction must be either bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND PRODUCT SPECIES MUST BE PLUS OR MINUS END
        product = products[1];
        if(product.find("PLUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
            }
            else {
                cout << "A plus end species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(product.find("MINUSEND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                ///get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
            }
            else {
                cout << "A minus species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Fifth species listed in a depolymerization reaction must be either plusend or minusend. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///Add depolymerization managers
        if(d == FilamentReactionDirection::FORWARD)
            _IFRxnManagers.emplace_back(new DepolyMinusEndManager(reactantTemplate, productTemplate, get<2>(r)));
        else
            _IFRxnManagers.emplace_back(new DepolyPlusEndManager(reactantTemplate, productTemplate, get<2>(r)));
    }

    
    ///set up reaction templates
    for(auto &r: chem.bindingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        ///read strings, and look up type
        
        ///Checks on number of reactants, products
        if(reactants.size() != 2 || products.size() != 1) {
            cout << "Invalid binding reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST REACTANT SPECIES MUST BE BOUND
        auto reactant = reactants[0];
        ///read strings, and look up type
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "First species listed in a binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND REACTANT SPECIES MUST BE BULK OR DIFFUSING
        reactant = reactants[1];
        if(reactant.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Second species listed in a binding reaction must be either bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        
        ///FIRST PRODUCT SPECIES MUST BE BOUND
        auto product = products[0];
        ///read strings, and look up type
        if(product.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Third species listed in a binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    
    
        ///add reaction
        _IFRxnManagers.emplace_back(new BasicBindingManager(reactantTemplate, productTemplate, get<2>(r)));
    
    }

    ///set up reaction templates
    for(auto &r: chem.unbindingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        ///read strings, and look up type
        
        ///Checks on number of reactants, products
        if(reactants.size() != 1 || products.size() != 2) {
            cout << "Invalid unbinding reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST REACTANT SPECIES MUST BE BOUND, LINKER, MOTOR
        auto reactant = reactants[0];
        ///read strings, and look up type
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(reactant.find("LINKER") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesLinker.begin(), _speciesLinker.end(), name);
            int position = 0;
            
            if(it != _speciesLinker.end()) {
                
                ///get position of iterator
                position = distance(_speciesLinker.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::LINKER));
            }
            else {
                cout << "A linker species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(reactant.find("MOTOR") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(it != _speciesMotor.end()) {
                
                ///get position of iterator
                position = distance(_speciesMotor.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
            }
            else {
                cout << "A motor species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "First species listed in an unbinding reaction must be bound, linker, or motor. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST PRODUCT SPECIES MUST BE BULK OR DIFFUSING
        auto product = products[0];
        if(product.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(product.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Second species listed in an unbinding reaction must be either bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        
        ///SECOND PRODUCT SPECIES MUST BE BOUND
        product = products[1];
        ///read strings, and look up type
        if(product.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Third species listed in an unbinding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        
        ///add reaction
        _IFRxnManagers.emplace_back(new UnbindingManager(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    for(auto &r: chem.motorWalkingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        ///read strings, and look up type
        ReactionType type;
        string species1;
        
        ///Checks on number of reactants, products
        if(reactants.size() != 2 || products.size() != 2) {
            cout << "Invalid motor walking reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST REACTANT SPECIES MUST BE MOTOR
        auto reactant = reactants[0];
        if(reactant.find("MOTOR") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(it != _speciesMotor.end()) {
                
                species1 = name;
                
                if(reactant.find("N+1") != string::npos)
                    type = ReactionType::MOTORWALKINGBACKWARD;
                else type = ReactionType::MOTORWALKINGFORWARD;
                
                ///get position of iterator
                position = distance(_speciesMotor.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
            }
            else {
                cout << "A motor species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "First species listed in a motor walking reaction must be motor. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        
        ///SECOND REACTANT SPECIES MUST BE BOUND
        reactant = reactants[1];
        ///read strings, and look up type
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Second species listed in a motor walking reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
   
        ///FIRST PRODUCT SPECIES MUST BE MOTOR
        auto product = products[0];
        if(product.find("MOTOR") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout << "Motor species in reactants and products of motor walking reaction must be same. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            if(it != _speciesMotor.end()) {
                
                ///get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
            }
            else {
                cout << "A motor species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Third species listed in a motor walking reaction must be motor. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///SECOND PRODUCT SPECIES MUST BE BOUND
        ///read strings, and look up type
        product = products[1];
        if(product.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout << "Fourth species listed in a motor walking reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        ///add reaction
        if(type == ReactionType::MOTORWALKINGFORWARD)
            _IFRxnManagers.emplace_back(new MotorWalkFManager(reactantTemplate, productTemplate, get<2>(r)));
        else
            _IFRxnManagers.emplace_back(new MotorWalkBManager(reactantTemplate, productTemplate, get<2>(r)));
    }
 
}

void SimpleManagerImpl::genCFRxnManagers(ChemistryData& chem) {
    
    int ID = 0;
    for(auto &r: chem.linkerBindingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        ReactionType type;
        string species1;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
    
        ///Checks on number of reactants, products
        if(reactants.size() != 3 || products.size() != 2) {
            cout << "Invalid linker binding reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST TWO SPECIES SHOULD BE BOUND
        auto reactant = reactants[0];
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "First species listed in a linker binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        reactant = reactants[1];
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Second species listed in a linker binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
        reactant = reactants[2];
        if(reactant.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Third species listed in a linker binding reaction must be bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        

        ///FIRST TWO SPECIES IN PRODUCTS MUST BE LINKER
        auto product = products[0];
        
        if(product.find("LINKER") != string::npos) {
            
            type = ReactionType::LINKERBINDING;
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesLinker.begin(), _speciesLinker.end(), name);
            int position = 0;
            
            species1 = name;
            
            if(it != _speciesLinker.end()) {
                
                ///get position of iterator
                position = distance(_speciesLinker.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::LINKER));
            }
            else {
                cout << "A linker species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Fourth species listed in a linker binding reaction must be linker. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        product = products[1];
        if(product.find("LINKER") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesLinker.begin(), _speciesLinker.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout << "Linker species in reactants and products of linker binding reaction must be same. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        
            if(it != _speciesLinker.end()) {
                
                ///get position of iterator
                position = distance(_speciesLinker.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::LINKER));
            }
            else {
                cout << "A linker species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else {
            cout << "Fifth species listed in a linker binding reaction must be linker. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    
        _CFRxnManagers.emplace_back(new LinkerBindingManager(reactantTemplate, productTemplate, get<2>(r), get<3>(r), get<4>(r), ID++));
    }
    
    for(auto &r: chem.motorBindingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        ReactionType type;
        string species1;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        ///Checks on number of reactants, products
        if(reactants.size() != 3 || products.size() != 2) {
            cout << "Invalid motor binding reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///FIRST TWO SPECIES SHOULD BE BOUND
        auto reactant = reactants[0];
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "First species listed in a motor binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        reactant = reactants[1];
        if(reactant.find("BOUND") != string::npos) {
            
            ///look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                ///get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
            }
            else {
                cout << "A bound species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Second species listed in a motor binding reaction must be bound. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        ///THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
        reactant = reactants[2];
        if(reactant.find("BULK") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                   [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            ///Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                   [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Third species listed in a cross filament binding reaction must be bulk or diffusing. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        
        ///FIRST TWO SPECIES IN PRODUCTS MUST BE MOTOR
        auto product = products[0];
        
        if(product.find("MOTOR") != string::npos) {
            
            type = ReactionType::MOTORBINDING;
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            species1 = name;
            
            if(it != _speciesMotor.end()) {
                
                ///get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
            }
            else {
                cout << "A motor species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "Fourth species listed in a motor binding reaction must be motor. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        product = products[1];
        if(product.find("MOTOR") != string::npos) {
            
            ///look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout << "Motor species in reactants and products of motor binding reaction must be same. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
 
            if(it != _speciesMotor.end()) {
                
                ///get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
            }
            else {
                cout << "A motor species that was included in a reaction was not initialized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else {
            cout << "Fifth species listed in a motor binding reaction must be motor. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        _CFRxnManagers.emplace_back(new MotorBindingManager(reactantTemplate, productTemplate, get<2>(r), get<3>(r), get<4>(r), ID++));
    }
}


void SimpleManagerImpl::genGeneralReactions(ChemistryData& chem, Compartment& protoCompartment) {
    
     ///go through reactions, add each
    for(auto &r: chem.genReactions) {
    
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                ///Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                       [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(CompartmentGrid::instance(compartmentGridKey())->findSpeciesBulkByName(name));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                ///Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                       [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
                if(it == chem.speciesDiffusing.end()) {
                    cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(protoCompartment.findSpeciesByName(name));
            }
            else {
                cout << "All reactants and products in a general reaction must be bulk or diffusing. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        for(auto &product : products) {
            if(product.find("BULK") != string::npos) {
                
                ///Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                       [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(CompartmentGrid::instance(compartmentGridKey())->findSpeciesBulkByName(name));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                ///Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                                       [name](tuple<string, int, double> element) { return get<0>(element) == name ? true : false; });
                if(it == chem.speciesDiffusing.end()) {
                    cout << "A diffusing species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(protoCompartment.findSpeciesByName(name));
            }
            else {
                cout << "All reactants and products in a general reaction must be bulk or diffusing. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        ///add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        
        ReactionBase* rxn;
        ///create the reaction
        
        ///<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r), true);
        ///<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r), true);
        ///<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r), true);
        ///<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r), true);
        ///<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r), true);
        ///<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r), true);
        ///<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r), true);
        ///<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r), true);
        else {
            cout << "General reaction specified does not match any existing templates. Exiting." <<endl;
            exit(EXIT_FAILURE);
        }
        
        ///add to compartment
        protoCompartment.addInternalReactionUnique(unique_ptr<ReactionBase>(rxn));
    }
}

void SimpleManagerImpl::genBulkReactions(ChemistryData& chem) {
    
    ///go through reactions, add each
    for(auto &r: chem.bulkReactions) {
        
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                ///Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                       [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(CompartmentGrid::instance(compartmentGridKey())->findSpeciesBulkByName(name));
            }
            else {
                cout << "All reactants and products in a bulk reaction must be bulk. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        for(auto &product : products) {
            if(product.find("BULK") != string::npos) {
                
                ///Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                       [name](tuple<string, int> element) { return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout << "A bulk species that was included in a reaction was not initialized. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(CompartmentGrid::instance(compartmentGridKey())->findSpeciesBulkByName(name));
            }
            else {
                cout << "All reactants and products in a bulk reaction must be bulk. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        }
        ///add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        ///create the reaction
        
        ///<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r));
        ///<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r));
        ///<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r));
        ///<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r));
        ///<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r));
        ///<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r));
        ///<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r));
        ///<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r));
        else {
            cout << "Bulk reaction specified does not match any existing templates. Exiting." <<endl;
            exit(EXIT_FAILURE);
        }
        
        ///add to grid
        CompartmentGrid::instance(compartmentGridKey())->addBulkReactionUnique(unique_ptr<ReactionBase>(rxn));
    }
}

void SimpleManagerImpl::initialize(ChemistryData& chem) {
    
    ///set static system ptr
    InternalFilamentRxnManager::_ps = _subSystem;
    CrossFilamentRxnManager::_ps = _subSystem;
    
    ///Copy all species from chem struct
    _speciesFilament =  chem.speciesFilament;
    _speciesPlusEnd  =  chem.speciesPlusEnd;
    _speciesMinusEnd =  chem.speciesMinusEnd;
    
    _speciesBound  =  chem.speciesBound;
    _speciesLinker =  chem.speciesLinker;
    _speciesMotor  =  chem.speciesMotor;
    
    ///Setup all species diffusing and bulk
    Compartment& cProto = CompartmentGrid::instance(compartmentGridKey())->getProtoCompartment();
    
    for(auto &sd : chem.speciesDiffusing)
        cProto.addSpeciesUnique(unique_ptr<Species>(
           new SpeciesDiffusing(get<0>(sd), get<1>(sd))), get<2>(sd));
    
    for(auto &sb : chem.speciesBulk)
        CompartmentGrid::instance(compartmentGridKey())->addSpeciesBulk(get<0>(sb), get<1>(sb));
    
    ///add reactions to protocompartment
    genGeneralReactions(chem, cProto);
    ///generate any bulk reactions
    genBulkReactions(chem);
    
    ///initialize all compartments equivalent to cproto
    for(auto &c : CompartmentGrid::instance(compartmentGridKey())->children()) {
        Compartment *C = static_cast<Compartment*>(c.get());
        *C = cProto;
    }
    for(auto &c : CompartmentGrid::instance(compartmentGridKey())->children()) {
        Compartment *C = static_cast<Compartment*>(c.get());
        C->generateAllDiffusionReactions();
    }
    
    ///add reactions to chemsim
    CompartmentGrid::instance(compartmentGridKey())->addChemSimReactions();
    
    ///create internal filament reaction managers
    genIFRxnManagers(chem);
    ///create cross filament reaction managers
    genCFRxnManagers(chem);
}

void SimpleManagerImpl::initializeCCylinder(CCylinder* cc, Filament *f,
                                            bool extensionFront, bool extensionBack, bool creation)
{
    ///maxlength is same length as mcylinder
    int maxlength = cc->getSize();
    Compartment* c = cc->getCompartment();
    
    ///add monomers to cylinder
    for(int i = 0; i < maxlength; i++) {
        
        CMonomer* m = new CMonomer();
        for(auto &f : _speciesFilament) {
            SpeciesFilament* sf =
                c->addSpeciesFilament(SpeciesNamesDB::Instance()->generateUniqueName(f));
            m->addSpeciesFilament(sf);
        }
        for (auto &p : _speciesPlusEnd) {
            SpeciesPlusEnd* sp =
                c->addSpeciesPlusEnd(SpeciesNamesDB::Instance()->generateUniqueName(p));
            m->addSpeciesPlusEnd(sp);
        }
        for (auto &mi : _speciesMinusEnd) {
            SpeciesMinusEnd* smi =
                c->addSpeciesMinusEnd(SpeciesNamesDB::Instance()->generateUniqueName(mi));
            m->addSpeciesMinusEnd(smi);
        }
        for (auto &b : _speciesBound) {
            SpeciesBound* sb =
                c->addSpeciesBound(SpeciesNamesDB::Instance()->generateUniqueName(b));
            m->addSpeciesBound(sb);
        }
        for (auto &l : _speciesLinker) {
            SpeciesLinker* sl =
                c->addSpeciesLinker(SpeciesNamesDB::Instance()->generateUniqueName(l));
            m->addSpeciesLinker(sl);
        }
        for (auto &mo : _speciesMotor) {
            SpeciesMotor* sm =
                c->addSpeciesMotor(SpeciesNamesDB::Instance()->generateUniqueName(mo));
            m->addSpeciesMotor(sm);
        }
        
        cc->addCMonomer(m);
    }
    ///get last ccylinder
    CCylinder* lastcc = nullptr;
 
    ///extension of front
    if(extensionFront) {
        lastcc = f->getCylinderVector().back()->getCCylinder();
        for(auto &r : _IFRxnManagers) r->addReaction(lastcc, cc);
    }
    ///extension of back
    else if(extensionBack) {
        lastcc = f->getCylinderVector().front()->getCCylinder();
        for(auto &r : _IFRxnManagers) r->addReaction(cc, lastcc);
    }

    else if(creation) {
        
        CMonomer* m2 = lastcc->getCMonomer(int(cc->getSize() / 2));
        CMonomer* m1 = lastcc->getCMonomer(int(cc->getSize() / 2) - 1);
        m2->speciesPlusEnd(0)->getRSpecies().setN(1);
        m1->speciesMinusEnd(0)->getRSpecies().setN(1);
    }
    
    ///Base case, initialization
    else {
        ///Check if this is the first cylinder
        if(f->getCylinderVector().size() != 0) {
            
            ///remove plus end from last, add to this.
            lastcc = f->getCylinderVector().back()->getCCylinder();
            CMonomer* m1 = lastcc->getCMonomer(lastcc->getSize() - 2);
            m1->speciesPlusEnd(0)->getRSpecies().setN(0);
            
            CMonomer* m2 = cc->getCMonomer(cc->getSize() - 2);
            m2->speciesPlusEnd(0)->getRSpecies().setN(1);
            m2->speciesBound(0)->getRSpecies().setN(1);
            
            ///fill last cylinder with default filament value
            for(int i = lastcc->getSize() - 2; i < lastcc->getSize(); i++) {
                lastcc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                lastcc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
                
            }
            ///fill new cylinder with default filament value
            for(int i = 0; i < cc->getSize() - 2; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
            for(auto &r : _IFRxnManagers) r->addReaction(lastcc, cc);
        }
        ///this is first one
        else {
            //set back and front
            CMonomer* m1 = cc->getCMonomer(cc->getSize() - 2);
            m1->speciesPlusEnd(0)->getRSpecies().setN(1);
            m1->speciesBound(0)->getRSpecies().setN(1);
            
            CMonomer* m2 = cc->getCMonomer(1);
            m2->speciesMinusEnd(0)->getRSpecies().setN(1);
            m2->speciesBound(0)->getRSpecies().setN(1);
            
            ///fill with default filament value
            for(int i = 2; i < cc->getSize() - 2; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
        }
    }    
    ///Add all reaction templates to this cylinder
    for(auto &r : _IFRxnManagers) { r->addReaction(cc); }

}

void SimpleManagerImpl::updateCCylinder(CCylinder* cc) {
    
    ///loop through all cross cylinder reactions, remove if cross filament
    auto ccReactions = cc->getCrossCylinderReactions();
    for(auto it = ccReactions.begin(); it != ccReactions.end(); it++) {
        
        auto ccOther = it->first;
        if(ccOther->getCylinder()->getFilament() != cc->getCylinder()->getFilament())
            cc->removeCrossCylinderReactions(ccOther);
    }
    
    ///Add reactions from manager, using the local neighbor list
    for(auto &r : _CFRxnManagers) {
        
        auto neighbors = r->getNeighborList()->getNeighbors(cc->getCylinder());
        for(auto &n : neighbors) r->addReaction(cc, static_cast<Cylinder*>(n)->getCCylinder());
    }
}


