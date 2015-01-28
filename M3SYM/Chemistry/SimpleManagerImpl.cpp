
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

#include <cmath>

#include "SimpleManagerImpl.h"

#include "ChemCallbacks.h"
#include "Parser.h"

#include "Bead.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

void SimpleManagerImpl::genIFRxnManagers(ChemistryData& chem) {
    
    //set up reaction templates
    for(auto &r: chem.polymerizationReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        FilamentReactionDirection d;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        //read strings, and look up type
        
        //Checks on number of reactants, products
        if(reactants.size() != POLYREACTANTS ||
           products.size() != POLYPRODUCTS - 1) {
            cout << "Invalid polymerization reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST SPECIES MUST BE BULK OR DIFFUSING
        auto reactant = reactants[0];
        if(reactant.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                              [name](tuple<string, int, string> element) {
                              return get<0>(element) == name ? true : false; });
                                       
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back( tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                              [name](tuple<string, int, double> element) {
                              return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout <<
            "First species listed in a polymerization reaction must be either bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND SPECIES MUST BE PLUS OR MINUS END
        reactant = reactants[1];
        if(reactant.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::PLUSEND));
                
                d = FilamentReactionDirection::FORWARD;
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(reactant.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::MINUSEND));
                
                d = FilamentReactionDirection::BACKWARD;
            }
            else {
                cout <<
                "A minus species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            
            cout <<
            "Second species listed in a polymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
            
        }
        
        //FIRST PRODUCT SPECIES MUST BE FILAMENT SPECIES
        auto product = products[0];
        if(product.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_speciesFilament.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                   SpeciesType::FILAMENT));
            }
            else {
                cout <<
                "A filament species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Third species listed in a polymerization reaction must be filament. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS OR MINUS END
        product = products[1];
        //read strings, and look up type
        if(product.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plusend species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(product.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                    SpeciesType::MINUSEND));
            }
            else {
                cout <<
                "A plusend species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Fourth species listed in a polymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //Add polymerization managers
        if(d == FilamentReactionDirection::FORWARD)
            _IFRxnManagers.emplace_back(
                new PolyPlusEndManager(reactantTemplate, productTemplate, get<2>(r)));
        else
            _IFRxnManagers.emplace_back(
                new PolyMinusEndManager(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    
    //set up reaction templates
    for(auto &r: chem.depolymerizationReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        FilamentReactionDirection d;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        //read strings, and look up type
        
        //Checks on number of reactants, products
        if(reactants.size() != DEPOLYREACTANTS - 1 ||
           products.size() != DEPOLYPRODUCTS) {
            cout << "Invalid depolymerization reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        //FIRST REACTANT SPECIES MUST BE FILAMENT SPECIES
        auto reactant = reactants[0];
        if(reactant.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_speciesFilament.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::FILAMENT));
            }
            else {
                cout <<
                "A filament species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "First species listed in a depolymerization reaction must be filament. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND REACTANT SPECIES MUST BE PLUS OR MINUS END
        reactant = reactants[1];
        //read strings, and look up type
        if(reactant.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::PLUSEND));
                
                d = FilamentReactionDirection::BACKWARD;
            }
            else {
                cout <<
                "A plusend species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(reactant.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::MINUSEND));
                d = FilamentReactionDirection::FORWARD;
            }
            else {
                cout <<
                "A minusend species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Second species listed in a depolymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        
        //FIRST PRODUCT SPECIES MUST BE BULK OR DIFFUSING
        auto product = products[0];
        if(product.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                              [name](tuple<string, int, string> element) {
                              return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(product.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                                [name](tuple<string, int, double> element) {
                                return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout <<
            "Third species listed in a depolymerization reaction must be either bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS OR MINUS END
        product = products[1];
        if(product.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                    SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else if(product.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                   SpeciesType::MINUSEND));
            }
            else {
                cout <<
                "A minus species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Fourth species listed in a depolymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //Add depolymerization managers
        if(d == FilamentReactionDirection::FORWARD)
            _IFRxnManagers.emplace_back(
            new DepolyMinusEndManager(reactantTemplate, productTemplate, get<2>(r)));
        else
            _IFRxnManagers.emplace_back(
            new DepolyPlusEndManager(reactantTemplate, productTemplate, get<2>(r)));
    }

    for(auto &r: chem.motorWalkingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        //read strings, and look up type
        ReactionType type;
        string species1;
        
        //Checks on number of reactants, products
        if(reactants.size() != MWALKINGREACTANTS ||
           products.size() != MWALKINGPRODUCTS) {
            cout << "Invalid motor walking reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST REACTANT SPECIES MUST BE MOTOR
        auto reactant = reactants[0];
        if(reactant.find("MOTOR") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(it != _speciesMotor.end()) {
                species1 = name;
                
                //check if forward or backward walking
                if(reactant.find("N+1") != string::npos)
                    type = ReactionType::MOTORWALKINGBACKWARD;
                else
                    type = ReactionType::MOTORWALKINGFORWARD;
                
                //get position of iterator
                position = distance(_speciesMotor.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                        SpeciesType::MOTOR));
            }
            else {
                cout <<
                "A motor species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "First species listed in a motor walking reaction must be motor. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }

        
        //SECOND REACTANT SPECIES MUST BE BOUND
        reactant = reactants[1];
        //read strings, and look up type
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                        SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "Second species listed in a motor walking reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
   
        //FIRST PRODUCT SPECIES MUST BE MOTOR
        auto product = products[0];
        if(product.find("MOTOR") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout <<
                "Motor species in reactants and products of motor walking reaction must be same. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            
            if((product.find("N+1") != string::npos && type != ReactionType::MOTORWALKINGFORWARD)) {
                
                cout <<
                "Motor walking reaction must have a direction (Check N and N+1 distinctions). Exiting."
                <<endl;
                exit(EXIT_FAILURE);
                
            }
            
            if(it != _speciesMotor.end()) {
                
                //get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                       SpeciesType::MOTOR));
            }
            else {
                cout <<
                "A motor species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Third species listed in a motor walking reaction must be motor. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND PRODUCT SPECIES MUST BE BOUND
        product = products[1];
        if(product.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                        SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "Fourth species listed in a motor walking reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }

        //add reaction
        if(type == ReactionType::MOTORWALKINGFORWARD)
            _IFRxnManagers.emplace_back(
                new MotorWalkFManager(reactantTemplate, productTemplate, get<2>(r)));
        else {
            _IFRxnManagers.emplace_back(
                new MotorWalkBManager(reactantTemplate, productTemplate, get<2>(r)));
        }
    }
    
    //set up reaction templates
    for(auto &r: chem.agingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        //read strings, and look up type
        
        //Checks on number of reactants, products
        if(reactants.size() != AGINGREACTANTS ||
           products.size() != AGINGPRODUCTS) {
            cout << "Invalid aging reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST REACTANT SPECIES MUST BE FILAMENT, PLUS OR MINUS END
        auto reactant = reactants[0];
        //read strings, and look up type
        if(reactant.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_speciesFilament.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::FILAMENT));
            }
            else {
                cout <<
                "A filament species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(reactant.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(reactant.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::MINUSEND));
            }
            else {
                cout <<
                "A minus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "First species listed in an aging reaction must be filament. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST PRODUCT SPECIES MUST BE FILAMENT, PLUS, OR MINUS END
        auto product = products[0];
        //read strings, and look up type
        if(product.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_speciesFilament.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                   SpeciesType::FILAMENT));
            }
            else {
                cout <<
                "A filament species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(product.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(product.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_speciesMinusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                    SpeciesType::MINUSEND));
            }
            else {
                cout <<
                "A minus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "Second species listed in an aging reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //add reaction
        _IFRxnManagers.emplace_back(
            new AgingManager(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    //set up reaction templates
    for(auto &r: chem.destructionReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        //read strings, and look up type
        
        //Checks on number of reactants, products
        if(reactants.size() != DESTRUCTIONREACTANTS ||
           products.size() != DESTRUCTIONPRODUCTS ) {
            cout << "Invalid destruction reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST SPECIES MUST BE PLUS END
        auto reactant = reactants[0];
        if(reactant.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            
            if(it != _speciesPlusEnd.end()) {
                //get position of iterator
                int position = distance(_speciesPlusEnd.begin(), it);
                reactantTemplate.push_back(
                    tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "First species listed in a destruction reaction must be plus end. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND SPECIES MUST BE MINUS END
        reactant = reactants[1];
        if(reactant.find("MINUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
            
            if(it != _speciesMinusEnd.end()) {
                //get position of iterator
                int position = distance(_speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(
                    tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
            }
            else {
                cout <<
                "A minus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Second species listed in a destruction reaction must be minus end. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }

        ///ALL PRODUCTS MUST BE BULK OR DIFFUSING
        for (auto &product : products) {
            
            if(product.find("BULK") != string::npos) {
                    
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                  [name](tuple<string, int, string> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                                  [name](tuple<string, int, double> element) {
                                  return get<0>(element) == name ? true : false; });
                if(it == chem.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
            }
            else {
                cout <<
                "Third species listed in a destruction reaction must be either bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
            
        //add reaction
        _IFRxnManagers.emplace_back(
            new DestructionManager(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    //set up reaction templates
    for(auto &r: chem.severingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        string reactant = get<0>(r);
        //read strings, and look up type
        
        
        // SPECIES MUST BE FILAMENT
        if(reactant.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
            int position = 0;
            
            if(it != _speciesFilament.end()) {
                //get position of iterator
                position = distance(_speciesFilament.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                    SpeciesType::FILAMENT));
            }
            else {
                cout <<
                "A filament species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Reactant species listed in a severing reaction must be filament. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //add reaction
        _IFRxnManagers.emplace_back(
          new SeveringManager(reactantTemplate, productTemplate, get<1>(r)));
    }
    
    //set up reaction templates
    for(auto &r: chem.branchingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        //Checks on number of reactants, products
        if(reactants.size() != BRANCHINGREACTANTS ||
           products.size() != BRANCHINGPRODUCTS) {
            cout << "Invalid branching reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST AND SECOND SPECIES MUST BE BULK OR DIFFUSING
        auto reactant = reactants[0];
        if(reactant.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                              [name](tuple<string, int, string> element) {
                              return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back( tuple<int, SpeciesType>(
            SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                              [name](tuple<string, int, double> element) {
                              return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
            SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout <<
            "First species listed in a branching reaction must be either bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        reactant = reactants[1];
        if(reactant.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                              [name](tuple<string, int, string> element) {
                                  return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back( tuple<int, SpeciesType>(
            SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                              [name](tuple<string, int, double> element) {
                                  return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
            SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout <<
            "Second species listed in a branching reaction must be either bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }

        //THIRD REACTANT SPECIES MUST BE BOUND
        reactant = reactants[2];
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "Third species listed in a branching reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }

        //FIRST PRODUCT SPECIES MUST BE BRANCHER
        auto product = products[0];
        if(product.find("BRANCHER") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesBrancher.begin(), _speciesBrancher.end(), name);
            int position = 0;
            
            if(it != _speciesBrancher.end()) {
                
                //get position of iterator
                position = distance(_speciesBrancher.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                   SpeciesType::BRANCHER));
            }
            else {
                cout <<
                "A brancher species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else{
            cout <<
            "Fourth species listed in a branching reaction must be brancher. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //SECOND PRODUCT SPECIES MUST BE PLUS END
        product = products[1];
        if(product.find("PLUSEND") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _speciesPlusEnd.end()) {
                //get position of iterator
                position = distance(_speciesPlusEnd.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                   SpeciesType::PLUSEND));
            }
            else {
                cout <<
                "A plus end species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Second product species listed in a branching reaction must be plus end. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //add reaction
        _IFRxnManagers.emplace_back(
        new BranchingManager(reactantTemplate, productTemplate, get<2>(r), get<3>(r)));
    }
}

void SimpleManagerImpl::genCFRxnManagers(ChemistryData& chem) {
    
    for(auto &r: chem.linkerReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        string species1;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
    
        //Checks on number of reactants, products
        if(reactants.size() != LMBINDINGREACTANTS ||
           products.size() != LMBINDINGPRODUCTS) {
            cout << "Invalid linker reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST TWO SPECIES SHOULD BE BOUND
        auto reactant = reactants[0];
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                       SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "First species listed in a linker reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        reactant = reactants[1];
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                        SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Second species listed in a linker reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
        reactant = reactants[2];
        if(reactant.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                [name](tuple<string, int, string> element) {
                                return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                              [name](tuple<string, int, double> element) {
                              return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout << "Third species listed in a linker reaction must be bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        

        //FIRST TWO SPECIES IN PRODUCTS MUST BE LINKER
        auto product = products[0];
        
        if(product.find("LINKER") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesLinker.begin(), _speciesLinker.end(), name);
            int position = 0;
            
            species1 = name;
            
            if(it != _speciesLinker.end()) {
                
                //get position of iterator
                position = distance(_speciesLinker.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::LINKER));
            }
            else {
                cout <<
                "A linker species that was included in a reaction was not initialized. Exiting." <<
                endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Fourth species listed in a linker reaction must be linker. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        product = products[1];
        if(product.find("LINKER") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesLinker.begin(), _speciesLinker.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout <<
                "Linker species in reactants and products of linker reaction must be same. Exiting." <<
                endl;
                exit(EXIT_FAILURE);
            }
        
            if(it != _speciesLinker.end()) {
                
                //get position of iterator
                position = distance(_speciesLinker.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::LINKER));
            }
            else {
                cout <<
                "A linker species that was included in a reaction was not initialized. Exiting." <<
                endl;
                exit(EXIT_FAILURE);
            }
        }
        
        else {
            cout <<
            "Fifth species listed in a linker reaction must be linker. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        //extend range of rMax, rMin by cylinder size
        double rMax = get<4>(r);
        double rMin = get<5>(r);
        _CFRxnManagers.emplace_back(
            new LinkerRxnManager(reactantTemplate, productTemplate,
                                 get<2>(r), get<3>(r), rMax, rMin));
    }
   
    for(auto &r: chem.motorReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        string species1;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        //Checks on number of reactants, products
        if(reactants.size() != LMBINDINGREACTANTS ||
           products.size() != LMBINDINGPRODUCTS) {
            cout << "Invalid motor reaction. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        //FIRST TWO SPECIES SHOULD BE BOUND
        auto reactant = reactants[0];
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                        SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "First species listed in a motor reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        reactant = reactants[1];
        if(reactant.find("BOUND") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_speciesBound.begin(), _speciesBound.end(), name);
            int position = 0;
            
            if(it != _speciesBound.end()) {
                
                //get position of iterator
                position = distance(_speciesBound.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                       SpeciesType::BOUND));
            }
            else {
                cout <<
                "A bound species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Second species listed in a motor reaction must be bound. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
        reactant = reactants[2];
        if(reactant.find("BULK") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                              [name](tuple<string, int, string> element) {
                              return get<0>(element) == name ? true : false; });
            
            if(it == chem.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(chem.speciesDiffusing.begin(),chem.speciesDiffusing.end(),
                              [name](tuple<string, int, double> element) {
                              return get<0>(element) == name ? true : false; });
            if(it == chem.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::Instance()->stringToInt(name), SpeciesType::DIFFUSING));
        }
        else {
            cout <<
            "Third species listed in a motor reaction must be bulk or diffusing. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        
        //FIRST TWO SPECIES IN PRODUCTS MUST BE MOTOR
        auto product = products[0];
        if(product.find("MOTOR") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            species1 = name;
            
            if(it != _speciesMotor.end()) {
                
                //get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                       SpeciesType::MOTOR));
            }
            else {
                cout <<
                "A motor species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Fourth species listed in a motor reaction must be motor. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        product = products[1];
        if(product.find("MOTOR") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_speciesMotor.begin(), _speciesMotor.end(), name);
            int position = 0;
            
            if(name != species1) {
                cout <<
                "Motor species in reactants and products of motor binding reaction must be same. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
 
            if(it != _speciesMotor.end()) {
                
                //get position of iterator
                position = distance(_speciesMotor.begin(), it);
                productTemplate.push_back(tuple<int, SpeciesType>(position,
                                                      SpeciesType::MOTOR));
            }
            else {
                cout <<
                "A motor species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout <<
            "Fifth species listed in a motor reaction must be motor. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        double rMax = get<4>(r);
        double rMin = get<5>(r);
        _CFRxnManagers.emplace_back(
            new MotorRxnManager(reactantTemplate, productTemplate,
                                get<2>(r), get<3>(r), rMax, rMin));
    }
}


void SimpleManagerImpl::copySpecies(ChemistryData& chem) {
    
    //Copy all species from chem struct, check lengths
    _speciesFilament =  chem.speciesFilament;
    if(_speciesFilament.size() < 1) {
        cout << "Must specify at least one filament species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    _speciesPlusEnd  =  chem.speciesPlusEnd;
    if(_speciesPlusEnd.size() < 1) {
        cout << "Must specify at least one plus end species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    if(_speciesPlusEnd.size() < _speciesFilament.size()) {
        cout << "There must be a plus end for every filament species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    
    _speciesMinusEnd =  chem.speciesMinusEnd;
    if(_speciesMinusEnd.size() < 1) {
        cout << "Must specify at least one minus end species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    if(_speciesMinusEnd.size() < _speciesFilament.size()) {
        cout << "There must be a minus end for every filament species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    
    _speciesBound  =  chem.speciesBound;
    if(_speciesBound.size() < 1) {
        cout << "Must specify at least one bound species. Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    _speciesLinker =  chem.speciesLinker;
    _speciesMotor  =  chem.speciesMotor;
    _speciesBrancher = chem.speciesBrancher;
}


void SimpleManagerImpl::genSpecies(ChemistryData& chem,
                                   Compartment& protoCompartment) {
    
    // add diffusing species (zero copy number for now)
    for(auto &sd : chem.speciesDiffusing)
        protoCompartment.addSpeciesUnique(unique_ptr<Species>(
        new SpeciesDiffusing(get<0>(sd), 0)), get<2>(sd));
    
    // add bulk species (with copy number)
    for(auto &sb : chem.speciesBulk)
        CompartmentGrid::instance()->
        addSpeciesBulk(get<0>(sb), get<1>(sb),
        (get<2>(sb) == "CONST") ? true : false);
}

void SimpleManagerImpl::initDiffusingCopyNumbers(ChemistryData& chem) {
    
    //look at copy number for each species
    for(auto &s : chem.speciesDiffusing) {
        
        auto name = get<0>(s);
        auto copyNumber = get<1>(s);
        
        //add randomly in compartment
        while (copyNumber > 0) {
            
            //find a random compartment within the boundary
            Compartment* randomCompartment;
            while(true) {
                randomCompartment = GController::getRandomCompartment();
                if(randomCompartment->isActivated()) break;
            }
            //find the species, increase copy number
            Species* species = randomCompartment->findSpeciesByName(name);
            species->up();
            copyNumber--;
        }
    }
}

void SimpleManagerImpl::genGeneralReactions(ChemistryData& chem,
                                            Compartment& protoCompartment) {
    
     //go through reactions, add each
    for(auto &r: chem.genReactions) {
    
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                 [name](tuple<string, int, string> element) {
                                 return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(
                CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it =
                    find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                            [name](tuple<string, int, double> element) {
                            return get<0>(element) == name ? true : false; });
                if(it == chem.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(protoCompartment.findSpeciesByName(name));
            }
            else {
                cout <<
                "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        for(auto &product : products) {
            if(product.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                  [name](tuple<string, int, string> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(
                CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it =
                    find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                            [name](tuple<string, int, double> element) {
                            return get<0>(element) == name ? true : false; });
                if(it == chem.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(protoCompartment.findSpeciesByName(name));
            }
            else {
                cout <<
                "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        //add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        //create the reaction
        
        //<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r), true);
        //<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r), true);
        //<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r), true);
        //<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r), true);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r), true);
        //<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r), true);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r), true);
        //<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r), true);
        else {
            cout <<
            "General reaction specified does not match any existing templates. Exiting."
            <<endl;
            exit(EXIT_FAILURE);
        }
        
        //add to compartment
        protoCompartment.addInternalReactionUnique(unique_ptr<ReactionBase>(rxn));
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void SimpleManagerImpl::genBulkReactions(ChemistryData& chem) {
    
    //go through reactions, add each
    for(auto &r: chem.bulkReactions) {
        
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                  [name](tuple<string, int, string> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(
                CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            else {
                cout <<
                "All reactants and products in a bulk reaction must be bulk. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        for(auto &product : products) {
            if(product.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                  [name](tuple<string, int, string> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == chem.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(
                CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            else {
                cout <<
                "All reactants and products in a bulk reaction must be bulk. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        //add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());
        
        ReactionBase* rxn;
        //create the reaction
        
        //<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r));
        //<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r));
        //<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r));
        //<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r));
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r));
        //<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r));
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r));
        //<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r));
        else {
            cout << "Bulk reaction specified does not match any existing templates. Exiting." <<endl;
            exit(EXIT_FAILURE);
        }
        
        //add to grid
        CompartmentGrid::instance()->
            addBulkReactionUnique(unique_ptr<ReactionBase>(rxn));
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void SimpleManagerImpl::genNucleationReactions(ChemistryData& chem) {
    
    //loop through all compartments
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = (Compartment*)(c.get());
        
        if(!C->isActivated()) continue;
        
        //go through reactions, add each
        for(auto &r: chem.nucleationReactions) {
            
            vector<Species*> reactantSpecies;
            
            vector<string> reactants = get<0>(r);
            vector<string> products = get<1>(r);
            
            if(reactants.size() != NUCLEATIONREACTANTS ||
               products.size() != NUCLEATIONPRODUCTS) {
                cout << "Invalid nucleation reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            bool diffusing = false;
            
            for(auto &reactant : reactants) {
                if(reactant.find("BULK") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chem.speciesBulk.begin(), chem.speciesBulk.end(),
                                      [name](tuple<string, int, string> element) {
                                      return get<0>(element) == name ? true : false; });
                    
                    if(it == chem.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(
                    CompartmentGrid::instance()->findSpeciesBulkByName(name));
                }
                
                else if(reactant.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it =
                    find_if(chem.speciesDiffusing.begin(), chem.speciesDiffusing.end(),
                            [name](tuple<string, int, double> element) {
                            return get<0>(element) == name ? true : false; });
                    if(it == chem.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));
                    diffusing = true;
                }
                else {
                    cout <<
                    "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            //add the reaction. The products will only be involved in creating the
            //callback needed to create a new filament
            ReactionBase* rxn = new Reaction<2,0>(reactantSpecies, get<2>(r));
            rxn->setReactionType(ReactionType::FILAMENTCREATION);
            
            C->addInternalReactionUnique(unique_ptr<ReactionBase>(rxn));
            
            reactantSpecies.clear();
            
            //now, loop through products, add callback
            short plusEnd;
            short minusEnd;
            short filament;
            
            //FIRST SPECIES MUST BE PLUS END
            auto product = products[0];
            if(product.find("PLUSEND") != string::npos) {
                
                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(_speciesPlusEnd.begin(), _speciesPlusEnd.end(), name);
                
                if(it != _speciesPlusEnd.end()) {
                    //get position of iterator
                    plusEnd = distance(_speciesPlusEnd.begin(), it);
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "First product species listed in a nucleation reaction must be plus end. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND SPECIES MUST BE FILAMENT
            product = products[1];
            if(product.find("FILAMENT") != string::npos) {
                
                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(_speciesFilament.begin(), _speciesFilament.end(), name);
                
                if(it != _speciesFilament.end()) {
                    //get position of iterator
                    filament = distance(_speciesFilament.begin(), it);
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Second product species listed in a nucleation reaction must be filament. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            
            //THIRD SPECIES MUST BE MINUSEND
            product = products[2];
            if(product.find("MINUSEND") != string::npos) {
                
                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(_speciesMinusEnd.begin(), _speciesMinusEnd.end(), name);
                
                if(it != _speciesMinusEnd.end()) {
                    //get position of iterator
                    minusEnd = distance(_speciesMinusEnd.begin(), it);
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Third product species listed in a nucleation reaction must be minus end. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            
            //if the reaction had any diffusing species, create the filament
            //in a random position within that compartment
            Compartment* creationCompartment;
            if(diffusing) creationCompartment = C;
            else creationCompartment = GController::getRandomCompartment();

            //now, add the callback
            FilamentCreationCallback
            fcallback(plusEnd, filament, minusEnd, _subSystem, creationCompartment);

            boost::signals2::shared_connection_block
            rcb(rxn->connect(fcallback,false));
        }
    }
}

void SimpleManagerImpl::initialize(ChemistryData& chem) {
    
    //set static system ptr
    InternalFilamentRxnManager::_ps = _subSystem;
    CrossFilamentRxnManager::_ps = _subSystem;
    
    //copy species
    copySpecies(chem);
    
    //Setup all species diffusing and bulk
    Compartment& cProto = CompartmentGrid::instance()->getProtoCompartment();
    
    //generate all species
    genSpecies(chem, cProto);
    
    //add reactions to protocompartment
    genGeneralReactions(chem, cProto);
    //generate any bulk reactions
    genBulkReactions(chem);
    
    //initialize all compartments equivalent to cproto
    //will copy all general and bulk reactions
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = (Compartment*)(c.get());
        *C = cProto;
    }
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = (Compartment*)(c.get());
        C->generateAllDiffusionReactions();
    }
    //generate nucleation reactions
    genNucleationReactions(chem);
    
    //intialize copy numbers of diffusing
    initDiffusingCopyNumbers(chem);
    
    //add reactions to chemsim
    CompartmentGrid::instance()->addChemSimReactions();
    
    //create internal filament reaction managers
    genIFRxnManagers(chem);
    //create cross filament reaction managers
    genCFRxnManagers(chem);
}

void SimpleManagerImpl::initializeCCylinder(CCylinder* cc, Filament *f,
                                            bool extensionFront,
                                            bool extensionBack,
                                            bool creation) {
    //maxlength is same length as mcylinder
    int maxlength = cc->getSize();
    Compartment* c = cc->getCompartment();
    
    //add monomers to cylinder
    for(int i = 0; i < maxlength; i++) {
        
        CMonomer* m = new CMonomer();
        for(auto &f : _speciesFilament) {
            SpeciesFilament* sf =
            c->addSpeciesFilament(
               SpeciesNamesDB::Instance()->genUniqueName(f));
            m->addSpeciesFilament(sf);
        }
        for (auto &p : _speciesPlusEnd) {
            SpeciesPlusEnd* sp =
            c->addSpeciesPlusEnd(
               SpeciesNamesDB::Instance()->genUniqueName(p));
            m->addSpeciesPlusEnd(sp);
        }
        for (auto &mi : _speciesMinusEnd) {
            SpeciesMinusEnd* smi =
            c->addSpeciesMinusEnd(
               SpeciesNamesDB::Instance()->genUniqueName(mi));
            m->addSpeciesMinusEnd(smi);
        }
        for (auto &b : _speciesBound) {
            SpeciesBound* sb =
            c->addSpeciesBound(
               SpeciesNamesDB::Instance()->genUniqueName(b));
            m->addSpeciesBound(sb);
        }
        for (auto &l : _speciesLinker) {
            SpeciesLinker* sl =
            c->addSpeciesLinker(
               SpeciesNamesDB::Instance()->genUniqueName(l));
            m->addSpeciesLinker(sl);
        }
        for (auto &mo : _speciesMotor) {
            SpeciesMotor* sm =
            c->addSpeciesMotor(
               SpeciesNamesDB::Instance()->genUniqueName(mo));
            m->addSpeciesMotor(sm);
        }
        for (auto &br : _speciesBrancher) {
            SpeciesBrancher* sb =
            c->addSpeciesBrancher(
               SpeciesNamesDB::Instance()->genUniqueName(br));
            m->addSpeciesBrancher(sb);
        }
        
        cc->addCMonomer(m);
    }
    //get last ccylinder
    CCylinder* lastcc = nullptr;
 
    //extension of front
    if(extensionFront) {
        lastcc = f->getCylinderVector().back()->getCCylinder();
        for(auto &r : _IFRxnManagers) r->addReaction(lastcc, cc);
    }
    //extension of back
    else if(extensionBack) {
        lastcc = f->getCylinderVector().front()->getCCylinder();
        for(auto &r : _IFRxnManagers) r->addReaction(cc, lastcc);
    }

    else if(creation) {
        //do nothing, this will be handled by the callback
    }
    
    //Base case, initialization
    else {
        //Check if this is the first cylinder
        if(f->getCylinderVector().size() != 0) {
            
            //remove plus end from last, add to this.
            lastcc = f->getCylinderVector().back()->getCCylinder();
            CMonomer* m1 = lastcc->getCMonomer(lastcc->getSize() - 1);
            m1->speciesPlusEnd(0)->down();
            
            CMonomer* m2 = cc->getCMonomer(cc->getSize() - 1);
            m2->speciesPlusEnd(0)->getRSpecies().setN(1);
            
            //fill last cylinder with default filament value
            m1->speciesFilament(0)->up();
            m1->speciesBound(0)->up();

            //fill new cylinder with default filament value
            for(int i = 0; i < cc->getSize() - 1; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
            for(auto &r : _IFRxnManagers) r->addReaction(lastcc, cc);
        }
        //this is first one
        else {
            //set back and front
            CMonomer* m1 = cc->getCMonomer(cc->getSize() - 1);
            m1->speciesPlusEnd(0)->getRSpecies().setN(1);
            
            CMonomer* m2 = cc->getCMonomer(0);
            m2->speciesMinusEnd(0)->getRSpecies().setN(1);
            
            //fill with default filament value
            for(int i = 1; i < cc->getSize() - 1; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->getRSpecies().setN(1);
                cc->getCMonomer(i)->speciesBound(0)->getRSpecies().setN(1);
            }
        }
    }    
    //Add all reaction templates to this cylinder
    for(auto &r : _IFRxnManagers) { r->addReaction(cc); }

}

void SimpleManagerImpl::updateCCylinder(CCylinder* cc) {
    
    //loop through all cross cylinder reactions, remove if cross filament
    auto ccReactions = cc->getCrossCylinderReactions();
    for(auto it = ccReactions.begin(); it != ccReactions.end(); it++) {
        
        auto ccOther = it->first;
        if(ccOther->getCylinder()->getFilament() != cc->getCylinder()->getFilament())
            cc->removeCrossCylinderReactions(ccOther, true);
    }
    
    //Add reactions from manager, using the local neighbor list
    for(auto &r : _CFRxnManagers) {
        
        auto neighbors = r->getNeighborList()->getNeighbors(cc->getCylinder());
        for(auto &n : neighbors)
            r->addReaction(cc, ((Cylinder*)(n))->getCCylinder());
    }
}


