
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include <cmath>

#include "SimpleManagerImpl.h"

#include "ChemCallbacks.h"
#include "Parser.h"

#include "Bead.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

void SimpleManagerImpl::initCMonomer(CMonomer* m, Compartment* c) {
    
    
    // FILAMENT SPECIES
    
    int fIndex = 0;
    for(auto &f : _chemData.speciesFilament) {
        SpeciesFilament* sf =
        c->addSpeciesFilament(
        SpeciesNamesDB::instance()->genUniqueName(f));
        m->_speciesFilament[fIndex] = sf;
        fIndex++;
    }
    for (auto &p : _chemData.speciesPlusEnd) {
        SpeciesPlusEnd* sp =
        c->addSpeciesPlusEnd(
        SpeciesNamesDB::instance()->genUniqueName(p));
        m->_speciesFilament[fIndex] = sp;
        fIndex++;
                            
    }
    for (auto &mi : _chemData.speciesMinusEnd) {
        SpeciesMinusEnd* smi =
        c->addSpeciesMinusEnd(
        SpeciesNamesDB::instance()->genUniqueName(mi));
        m->_speciesFilament[fIndex] = smi;
        fIndex++;
    }
                            
    // BOUND SPECIES
    
    int bIndex = 0;
    for (auto &b : _chemData.speciesBound) {
        SpeciesBound* sb =
        c->addSpeciesBound(
        SpeciesNamesDB::instance()->genUniqueName(b));
        m->_speciesBound[bIndex] = sb;
        bIndex++;
    }
    for (auto &l : _chemData.speciesLinker) {
        SpeciesLinker* sl =
        c->addSpeciesLinker(
        SpeciesNamesDB::instance()->genUniqueName(l));
        m->_speciesBound[bIndex] = sl;
        bIndex++;
    }
    for (auto &mo : _chemData.speciesMotor) {
        SpeciesMotor* sm =
        c->addSpeciesMotor(
        SpeciesNamesDB::instance()->genUniqueName(mo));
        m->_speciesBound[bIndex] = sm;
        bIndex++;
    }
    for (auto &br : _chemData.speciesBrancher) {
        SpeciesBrancher* sbr =
        c->addSpeciesBrancher(
        SpeciesNamesDB::instance()->genUniqueName(br));
        m->_speciesBound[bIndex] = sbr;
        bIndex++;
    }
}

void SimpleManagerImpl::genFilRxnTemplates() {
    
    //set up reaction templates
    for(auto &r: _chemData.polymerizationReactions) {
        
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
            auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                              [name](tuple<string, int, string, double> element) {
                              return get<0>(element) == name ? true : false; });
                                       
            if(it == _chemData.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back( tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(reactant.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                              [name](tuple<string, int, double, double> element) {
                              return get<0>(element) == name ? true : false; });
            if(it == _chemData.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::DIFFUSING));
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
                reactantTemplate.push_back(tuple<int, SpeciesType>(position,
                                                     SpeciesType::MINUSEND));
                
                d = FilamentReactionDirection::BACKWARD;
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
            "Second species listed in a polymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
            
        }
        
        //FIRST PRODUCT SPECIES MUST BE FILAMENT SPECIES
        auto product = products[0];
        if(product.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesFilament.begin(), it);
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
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
            _filRxnTemplates.emplace_back(
            new PolyPlusEndTemplate(reactantTemplate, productTemplate, get<2>(r)));
        else
            _filRxnTemplates.emplace_back(
            new PolyMinusEndTemplate(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    
    //set up reaction templates
    for(auto &r: _chemData.depolymerizationReactions) {
        
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
            auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesFilament.begin(), it);
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
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
            auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                              [name](tuple<string, int, string, double> element) {
                              return get<0>(element) == name ? true : false; });
            
            if(it == _chemData.speciesBulk.end()) {
                cout <<
                "A bulk species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::BULK));
        }
        
        else if(product.find("DIFFUSING") != string::npos) {
            
            //Look up species, make sure in list
            string name = product.substr(0, product.find(":"));
            auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                [name](tuple<string, int, double, double> element) {
                                return get<0>(element) == name ? true : false; });
            if(it == _chemData.speciesDiffusing.end()) {
                cout <<
                "A diffusing species that was included in a reaction was not initialized. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::DIFFUSING));
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
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
        else {
            cout <<
            "Fourth species listed in a depolymerization reaction must be either plusend or minusend. Exiting."
            << endl;
            exit(EXIT_FAILURE);
        }
        
        //Add depolymerization managers
        if(d == FilamentReactionDirection::FORWARD)
            _filRxnTemplates.emplace_back(
            new DepolyMinusEndTemplate(reactantTemplate, productTemplate, get<2>(r)));
        else
            _filRxnTemplates.emplace_back(
            new DepolyPlusEndTemplate(reactantTemplate, productTemplate, get<2>(r)));
    }

    for(auto &r: _chemData.motorWalkingReactions) {
        
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
            auto it = find(_chemData.speciesMotor.begin(), _chemData.speciesMotor.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMotor.end()) {
                species1 = name;
                
                //check if forward or backward walking
                if(reactant.find("N+1") != string::npos)
                    type = ReactionType::MOTORWALKINGBACKWARD;
                else
                    type = ReactionType::MOTORWALKINGFORWARD;
                
                //get position of iterator
                position = distance(_chemData.speciesMotor.begin(), it);
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
            auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesBound.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesBound.begin(), it);
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
            auto it = find(_chemData.speciesMotor.begin(), _chemData.speciesMotor.end(), name);
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
            
            if(it != _chemData.speciesMotor.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMotor.begin(), it);
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
            auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesBound.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesBound.begin(), it);
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
            _filRxnTemplates.emplace_back(
            new MotorWalkFTemplate(reactantTemplate, productTemplate, get<2>(r)));
        else {
            _filRxnTemplates.emplace_back(
            new MotorWalkBTemplate(reactantTemplate, productTemplate, get<2>(r)));
        }
    }
    
    //set up reaction templates
    for(auto &r: _chemData.agingReactions) {
        
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
            auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesFilament.begin(), it);
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
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
            auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesFilament.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesFilament.begin(), it);
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesPlusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesMinusEnd.end()) {
                
                //get position of iterator
                position = distance(_chemData.speciesMinusEnd.begin(), it);
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
        _filRxnTemplates.emplace_back(
        new AgingTemplate(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    
    //set up reaction templates
    for(auto &r: _chemData.destructionReactions) {
        
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
            auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
            
            if(it != _chemData.speciesPlusEnd.end()) {
                //get position of iterator
                int position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
            
            if(it != _chemData.speciesMinusEnd.end()) {
                //get position of iterator
                int position = distance(_chemData.speciesMinusEnd.begin(), it);
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double> element) {
                                  return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::DIFFUSING));
            }
            else {
                cout <<
                "Third species listed in a destruction reaction must be either bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
            
        //add reaction
        _filRxnTemplates.emplace_back(
        new DestructionTemplate(reactantTemplate, productTemplate, get<2>(r)));
    }
    
    //set up reaction templates
    for(auto &r: _chemData.severingReactions) {
        
        vector<tuple<int, SpeciesType>> reactantTemplate;
        vector<tuple<int, SpeciesType>> productTemplate;
        
        string reactant = get<0>(r);
        //read strings, and look up type
        
        
        // SPECIES MUST BE FILAMENT
        if(reactant.find("FILAMENT") != string::npos) {
            
            //look up species, make sure in list
            string name = reactant.substr(0, reactant.find(":"));
            auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
            int position = 0;
            
            if(it != _chemData.speciesFilament.end()) {
                //get position of iterator
                position = distance(_chemData.speciesFilament.begin(), it);
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
        _filRxnTemplates.emplace_back(
        new SeveringTemplate(reactantTemplate, productTemplate, get<1>(r)));
    }
}

void SimpleManagerImpl::genFilBindingManagers() {
    
    //loop through all compartments
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = (Compartment*)(c.get());
        
        if(!C->isActivated()) continue;
    
        for(auto &r: _chemData.branchingReactions) {
            
            vector<Species*> reactantSpecies;
            vector<Species*> productSpecies;
            
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double> element) {
                                      return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(C->findSpeciesByName(name));
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(CompartmentGrid::instance()->findSpeciesBulkByName(name));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double> element) {
                                      return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(C->findSpeciesByName(name));
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
                auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound.begin(), it);
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
                auto it = find(_chemData.speciesBrancher.begin(), _chemData.speciesBrancher.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBrancher.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBrancher.begin(), it);
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
                auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd.end()) {
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd.begin(), it);
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
            
            //Create reaction
            .emplace_back(
            new BranchingManager(reactantTemplate, productTemplate, get<2>(r), get<3>(r)));
        }
        
        
        for(auto &r: _chemData.linkerReactions) {
            
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
                auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound.begin(), it);
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
                auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound.begin(), it);
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                    [name](tuple<string, int, string, double> element) {
                                    return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double> element) {
                                  return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::DIFFUSING));
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
                auto it = find(_chemData.speciesLinker.begin(), _chemData.speciesLinker.end(), name);
                int position = 0;
                
                species1 = name;
                
                if(it != _chemData.speciesLinker.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesLinker.begin(), it);
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
                auto it = find(_chemData.speciesLinker.begin(), _chemData.speciesLinker.end(), name);
                int position = 0;
                
                if(name != species1) {
                    cout <<
                    "Linker species in reactants and products of linker reaction must be same. Exiting." <<
                    endl;
                    exit(EXIT_FAILURE);
                }
            
                if(it != _chemData.speciesLinker.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesLinker.begin(), it);
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
       
        for(auto &r: _chemData.motorReactions) {
            
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
                auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound.begin(), it);
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
                auto it = find(_chemData.speciesBound.begin(), _chemData.speciesBound.end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound.begin(), it);
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::BULK));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double> element) {
                                  return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::instance()->stringToInt(name), SpeciesType::DIFFUSING));
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
                auto it = find(_chemData.speciesMotor.begin(), _chemData.speciesMotor.end(), name);
                int position = 0;
                
                species1 = name;
                
                if(it != _chemData.speciesMotor.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMotor.begin(), it);
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
                auto it = find(_chemData.speciesMotor.begin(), _chemData.speciesMotor.end(), name);
                int position = 0;
                
                if(name != species1) {
                    cout <<
                    "Motor species in reactants and products of motor binding reaction must be same. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
     
                if(it != _chemData.speciesMotor.end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMotor.begin(), it);
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
}


void SimpleManagerImpl::genSpecies(Compartment& protoCompartment) {
    
    // add diffusing species (zero copy number for now)
    for(auto &sd : _chemData.speciesDiffusing)
        protoCompartment.addSpeciesUnique(unique_ptr<Species>(
        new SpeciesDiffusing(get<0>(sd), 0)), get<2>(sd));
    
    // add bulk species (zero copy number for now)
    for(auto &sb : _chemData.speciesBulk)
        CompartmentGrid::instance()->
        addSpeciesBulk(get<0>(sb), 0, (get<2>(sb) == "CONST") ? true : false);
}

void SimpleManagerImpl::updateCopyNumbers() {

    //look at copy number for each species
    for(auto &s : _chemData.speciesDiffusing) {
        
        auto name = get<0>(s);
        auto copyNumber = get<1>(s);
        auto releaseTime = get<3>(s);
        
        if(tau() >= releaseTime) {
        
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
            
            //set zero copy number
            get<1>(s) = 0;
        }
    }
    
    for(auto &s : _chemData.speciesBulk) {
        
        auto name = get<0>(s);
        auto copyNumber = get<1>(s);
        auto releaseTime = get<3>(s);
        
        if(tau() >= releaseTime && copyNumber != 0) {
        
            //find the species, set copy number
            Species* species = CompartmentGrid::instance()->
                               findSpeciesBulkByName(name);
            species->setN(copyNumber);
            
            //update reactions
            species->getRSpecies().activateAssocReactantReactions();
            
            //set zero copy number
            get<1>(s) = 0;
        }
    }
}

void SimpleManagerImpl::genGeneralReactions(Compartment& protoCompartment) {
    
     //go through reactions, add each
    for(auto &r: _chemData.genReactions) {
    
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                 [name](tuple<string, int, string, double> element) {
                                 return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
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
                    find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                            [name](tuple<string, int, double, double> element) {
                            return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
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
                    find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                            [name](tuple<string, int, double, double> element) {
                            return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
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

void SimpleManagerImpl::genBulkReactions() {
    
    //go through reactions, add each
    for(auto &r: _chemData.bulkReactions) {
        
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
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
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, string, double> element) {
                                  return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
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
            cout << "Bulk reaction specified does not match any existing templates. Exiting."
            <<endl;
            exit(EXIT_FAILURE);
        }
        
        //add to grid
        CompartmentGrid::instance()->
            addBulkReactionUnique(unique_ptr<ReactionBase>(rxn));
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void SimpleManagerImpl::genNucleationReactions() {
    
#if !defined(REACTION_SIGNALING)
    if(!_chemData.nucleationReactions.empty()) {
        
        cout << "Nucleation reactions rely on reaction signaling. Please set this "
             << "compilation macro. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
#endif
    
    //loop through all compartments
    for(auto &c : CompartmentGrid::instance()->children()) {
        Compartment *C = (Compartment*)(c.get());
        
        if(!C->isActivated()) continue;
        
        //go through reactions, add each
        for(auto &r: _chemData.nucleationReactions) {
            
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
                    auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                      [name](tuple<string, int, string, double> element) {
                                      return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
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
                    find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                            [name](tuple<string, int, double, double> element) {
                            return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
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
                auto it = find(_chemData.speciesPlusEnd.begin(), _chemData.speciesPlusEnd.end(), name);
                
                if(it != _chemData.speciesPlusEnd.end()) {
                    //get position of iterator
                    plusEnd = distance(_chemData.speciesPlusEnd.begin(), it);
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
                auto it = find(_chemData.speciesFilament.begin(), _chemData.speciesFilament.end(), name);
                
                if(it != _chemData.speciesFilament.end()) {
                    //get position of iterator
                    filament = distance(_chemData.speciesFilament.begin(), it);
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
                auto it = find(_chemData.speciesMinusEnd.begin(), _chemData.speciesMinusEnd.end(), name);
                
                if(it != _chemData.speciesMinusEnd.end()) {
                    //get position of iterator
                    minusEnd = distance(_chemData.speciesMinusEnd.begin(), it);
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
#ifdef REACTION_SIGNALING
            FilamentCreationCallback
            fcallback(plusEnd, filament, minusEnd, _subSystem, creationCompartment);

            boost::signals2::shared_connection_block
            rcb(rxn->connect(fcallback,false));
#endif
        }
    }
}


void SimpleManagerImpl::configureCMonomer() {
    
    //set up static CMonomer things
    CMonomer::_numFSpecies = _chemData.speciesFilament.size() +
                             _chemData.speciesPlusEnd.size()  +
                             _chemData.speciesMinusEnd.size();
    
    CMonomer::_numBSpecies = _chemData.speciesBound.size()   +
                             _chemData.speciesLinker.size()  +
                             _chemData.speciesMotor.size()   +
                             _chemData.speciesBrancher.size();
    
    //set up species offsets
    short o1 = _chemData.speciesFilament.size();
    short o2 = o1 + _chemData.speciesPlusEnd.size();
    
    short o3 = _chemData.speciesBound.size();
    short o4 = o3 + _chemData.speciesLinker.size();
    short o5 = o4 + _chemData.speciesMotor.size();
    
    //create offset vector for filament
    CMonomer::_speciesFilamentIndex.push_back(0);
    CMonomer::_speciesFilamentIndex.push_back(o1);
    CMonomer::_speciesFilamentIndex.push_back(o2);
    
    //create offset vector for bound
    CMonomer::_speciesBoundIndex.push_back(0);
    CMonomer::_speciesBoundIndex.push_back(o3);
    CMonomer::_speciesBoundIndex.push_back(o4);
    CMonomer::_speciesBoundIndex.push_back(o5);
}


void SimpleManagerImpl::initializeSystem() {
    
    //set static system ptr
    FilamentReactionTemplate::_ps = _subSystem;
    
    //init rand number gen
    FilamentBindingManager::_eng =
    new mt19937(static_cast<unsigned long>(time(nullptr)));
    
    //config CMonomer
    configureCMonomer();
    
    //Setup all species diffusing and bulk
    Compartment& cProto = CompartmentGrid::instance()->
                          getProtoCompartment();
    
    genSpecies(cProto);
    
    //will print reactions as well
    genGeneralReactions(cProto);
    genBulkReactions();
    
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
    //try initial copy number setting
    updateCopyNumbers();
    
    //generate the nucleation reactions in system
    genNucleationReactions();
    //generate filament binding managers
    genFilBindingManagers();

    //add reactions to chemsim
    CompartmentGrid::instance()->addChemSimReactions();
    
    //create filament reaction templates
    genFilRxnTemplates();
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
        initCMonomer(m, c);
        cc->addCMonomer(m);
    }
    
    //get last ccylinder
    CCylinder* lastcc = nullptr;
 
    //extension of front
    if(extensionFront) {
        lastcc = f->getCylinderVector().back()->getCCylinder();
        for(auto &r : _filRxnTemplates) r->addReaction(lastcc, cc);
    }
    //extension of back
    else if(extensionBack) {
        lastcc = f->getCylinderVector().front()->getCCylinder();
        for(auto &r : _filRxnTemplates) r->addReaction(cc, lastcc);
    }

    else if(creation) {
        //do nothing, this will be handled by the callback
    }
    
    //Base case, initialization
    else {
        //Check if this is the first cylinder
        if(!f->getCylinderVector().empty()) {
            
            //remove plus end from last, add to this.
            lastcc = f->getCylinderVector().back()->getCCylinder();
            CMonomer* m1 = lastcc->getCMonomer(lastcc->getSize() - 1);
            m1->speciesPlusEnd(0)->down();
            
            CMonomer* m2 = cc->getCMonomer(cc->getSize() - 1);
            m2->speciesPlusEnd(0)->up();
            
            //fill last cylinder with default filament value
            m1->speciesFilament(0)->up();
            m1->speciesBound(BOUND_EMPTY)->up();

            //fill new cylinder with default filament value
            for(int i = 0; i < cc->getSize() - 1; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->up();
                cc->getCMonomer(i)->speciesBound(BOUND_EMPTY)->up();
            }
            for(auto &r : _filRxnTemplates) r->addReaction(lastcc, cc);
        }
        //this is first one
        else {
            //set back and front
            CMonomer* m1 = cc->getCMonomer(cc->getSize() - 1);
            m1->speciesPlusEnd(0)->up();
            
            CMonomer* m2 = cc->getCMonomer(0);
            m2->speciesMinusEnd(0)->up();
            
            //fill with default filament value
            for(int i = 1; i < cc->getSize() - 1; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->up();
                cc->getCMonomer(i)->speciesBound(BOUND_EMPTY)->up();
            }
        }
    }    
    //Add all reaction templates to this cylinder
    for(auto &r : _filRxnTemplates) { r->addReaction(cc); }
}

