
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

#include <cmath>

#include "ChemManager.h"

#include "ChemCallbacks.h"
#include "CompartmentGrid.h"
#include "ReactionTemplate.h"
#include "BindingManager.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

void ChemManager::setupBindingSites() {
    
    //set binding indices
    //check if binding sites are valid and mark
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        
        if(_chemData.B_BINDING_INDEX[filType] != "") {
            auto it = find(_chemData.speciesBound[filType].begin(),
                           _chemData.speciesBound[filType].end(),
                           _chemData.B_BINDING_INDEX[filType]);
            
            if(it == _chemData.speciesBound[filType].end()) {
                
                cout << "The brancher binding site listed is not a valid bound species. Exiting."
                     << endl;
                 exit(EXIT_FAILURE);
            }
            else {
                SysParams::CParams.brancherBoundIndex[filType] = it - _chemData.speciesBound[filType].begin();
            }
        }
        
        if(_chemData.L_BINDING_INDEX[filType] != "") {
            
            auto it = find(_chemData.speciesBound[filType].begin(),
                           _chemData.speciesBound[filType].end(),
                           _chemData.L_BINDING_INDEX[filType]);
            
            if(it == _chemData.speciesBound[filType].end()) {
                
                cout << "The linker binding site listed is not a valid bound species. Exiting."
                     << endl;
                 exit(EXIT_FAILURE);
            }
            else {
                SysParams::CParams.linkerBoundIndex[filType] = it - _chemData.speciesBound[filType].begin();
            }
        }
        
        if(_chemData.M_BINDING_INDEX[filType] != "") {
            
            auto it = find(_chemData.speciesBound[filType].begin(),
                           _chemData.speciesBound[filType].end(),
                           _chemData.M_BINDING_INDEX[filType]);
            
            if(it == _chemData.speciesBound[filType].end()) {
                
                cout << "The motor binding site listed is not a valid bound species. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else {
                SysParams::CParams.motorBoundIndex[filType] = it - _chemData.speciesBound[filType].begin();
            }
        }
        
        //for initialization of cylinders
        SysParams::CParams.bindingIndices[filType].push_back(SysParams::CParams.brancherBoundIndex[filType]);
        
        if(SysParams::CParams.brancherBoundIndex[filType] != SysParams::CParams.linkerBoundIndex[filType])
             SysParams::CParams.bindingIndices[filType].push_back(SysParams::CParams.linkerBoundIndex[filType]);
        
        if(SysParams::CParams.brancherBoundIndex[filType] != SysParams::CParams.motorBoundIndex[filType] &&
           SysParams::CParams.linkerBoundIndex[filType]   != SysParams::CParams.motorBoundIndex[filType])
            SysParams::CParams.bindingIndices[filType].push_back(SysParams::CParams.motorBoundIndex[filType]);
    }
}

void ChemManager::configCMonomer() {
    
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        //set up static CMonomer things
        CMonomer::_numFSpecies[filType] = _chemData.speciesFilament[filType].size() +
                                        _chemData.speciesPlusEnd[filType].size()  +
                                        _chemData.speciesMinusEnd[filType].size();
        
        CMonomer::_numBSpecies[filType] = _chemData.speciesBound[filType].size()   +
                                        _chemData.speciesLinker[filType].size()  +
                                        _chemData.speciesMotor[filType].size()   +
                                        _chemData.speciesBrancher[filType].size();
        
        //set up species offsets
        short o1 = _chemData.speciesFilament[filType].size();
        short o2 = o1 + _chemData.speciesPlusEnd[filType].size();
        
        short o3 = _chemData.speciesBound[filType].size();
        short o4 = o3 + _chemData.speciesLinker[filType].size();
        short o5 = o4 + _chemData.speciesMotor[filType].size();
        
        //create offset vector for filament
        CMonomer::_speciesFilamentIndex[filType].insert(
        CMonomer::_speciesFilamentIndex[filType].end(), {0,o1,o2});
        
        //create offset vector for bound
        CMonomer::_speciesBoundIndex[filType].insert(
        CMonomer::_speciesBoundIndex[filType].end(), {0,o3,o4,o5});
    }
}

void ChemManager::initCMonomer(CMonomer* m, short filamentType, Compartment* c) {
    
    // FILAMENT SPECIES
    
    int fIndex = 0;
    for(auto &f : _chemData.speciesFilament[filamentType]) {
        SpeciesFilament* sf =
        c->addSpeciesFilament(SpeciesNamesDB::genUniqueFilName(f));
        m->_speciesFilament[fIndex] = sf;
        fIndex++;
    }
    for (auto &p : _chemData.speciesPlusEnd[filamentType]) {
        SpeciesPlusEnd* sp =
        c->addSpeciesPlusEnd(SpeciesNamesDB::genUniqueFilName(p));
        m->_speciesFilament[fIndex] = sp;
        fIndex++;
        
    }
    for (auto &mi : _chemData.speciesMinusEnd[filamentType]) {
        SpeciesMinusEnd* smi =
        c->addSpeciesMinusEnd(SpeciesNamesDB::genUniqueFilName(mi));
        m->_speciesFilament[fIndex] = smi;
        fIndex++;
    }
    
    // BOUND SPECIES
    
    int bIndex = 0;
    for (auto &b : _chemData.speciesBound[filamentType]) {
        SpeciesBound* sb =
        c->addSpeciesBound(SpeciesNamesDB::genUniqueFilName(b));
        m->_speciesBound[bIndex] = sb;
        bIndex++;
    }
    for (auto &l : _chemData.speciesLinker[filamentType]) {
        SpeciesLinker* sl =
        c->addSpeciesLinker(SpeciesNamesDB::genUniqueFilName(l));
        m->_speciesBound[bIndex] = sl;
        bIndex++;
    }
    for (auto &mo : _chemData.speciesMotor[filamentType]) {
        SpeciesMotor* sm =
        c->addSpeciesMotor(SpeciesNamesDB::genUniqueFilName(mo));
        m->_speciesBound[bIndex] = sm;
        bIndex++;
    }
    for (auto &br : _chemData.speciesBrancher[filamentType]) {
        SpeciesBrancher* sbr =
        c->addSpeciesBrancher(SpeciesNamesDB::genUniqueFilName(br));
        m->_speciesBound[bIndex] = sbr;
        bIndex++;
    }
}

void ChemManager::genFilReactionTemplates() {
    
    
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        //set up reaction templates
        for(auto &r: _chemData.polymerizationReactions[filType]) {
            
            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;
            FilamentReactionDirection d;
            
            vector<string> reactants = get<0>(r);
            vector<string> products = get<1>(r);
            //read strings, and look up type
            
            //Checks on number of reactants, products
            if(reactants.size() != POLYREACTANTS ||
               products.size() != POLYPRODUCTS) {
                cout << "Invalid polymerization reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            //FIRST SPECIES MUST BE BULK OR DIFFUSING
            auto reactant = reactants[0];
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, double, double, string> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back( tuple<int, SpeciesType>(
                SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double, double, string, int> element) {
                                      return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                    
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                    
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
                auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesFilament[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesFilament[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
                _filRxnTemplates[filType].emplace_back(
                new PolyPlusEndTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
            else
                _filRxnTemplates[filType].emplace_back(
                new PolyMinusEndTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
        }
        
        //set up reaction templates
        for(auto &r: _chemData.depolymerizationReactions[filType]) {
            
            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;
            FilamentReactionDirection d;
            
            vector<string> reactants = get<0>(r);
            vector<string> products = get<1>(r);
            //read strings, and look up type
            
            //Checks on number of reactants, products
            if(reactants.size() != DEPOLYREACTANTS ||
               products.size() != DEPOLYPRODUCTS) {
                cout << "Invalid depolymerization reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            //FIRST REACTANT SPECIES MUST BE FILAMENT SPECIES
            auto reactant = reactants[0];
            if(reactant.find("FILAMENT") != string::npos) {
                
                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesFilament[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                    
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
                                  [name](tuple<string, int, double, double, string> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                  [name](tuple<string, int, double, double, double, string, int> element) {
                                      return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
                _filRxnTemplates[filType].emplace_back(
                new DepolyMinusEndTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
            else
                _filRxnTemplates[filType].emplace_back(
                new DepolyPlusEndTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
        }
        
        for(auto &r: _chemData.motorWalkingReactions[filType]) {
            
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
                auto it = find(_chemData.speciesMotor[filType].begin(), _chemData.speciesMotor[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMotor[filType].end()) {
                    species1 = name;
                    
                    //check if forward or backward walking
                    if(reactant.find("N+1") != string::npos)
                        type = ReactionType::MOTORWALKINGBACKWARD;
                    else
                        type = ReactionType::MOTORWALKINGFORWARD;
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMotor[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
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
            
            
            //SECOND REACTANT SPECIES MUST BE EMPTY SITE
            reactant = reactants[1];
            //read strings, and look up type
            if(reactant.find("BOUND") != string::npos) {
                
                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound[filType].begin(), it);
                    
                    if(position != SysParams::CParams.motorBoundIndex[filType]) {
                        cout <<
                        "Second species listed in a motor walking reaction must be the corresponding motor empty site. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
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
                auto it = find(_chemData.speciesMotor[filType].begin(), _chemData.speciesMotor[filType].end(), name);
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
                
                if(it != _chemData.speciesMotor[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMotor[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
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
            
            //SECOND PRODUCT SPECIES MUST BE EMPTY SITE
            product = products[1];
            if(product.find("BOUND") != string::npos) {
                
                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesBound[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesBound[filType].begin(), it);
                    
                    if(position != SysParams::CParams.motorBoundIndex[filType]) {
                        cout <<
                        "Second species listed in a motor walking reaction must be the corresponding motor empty site. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
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
                _filRxnTemplates[filType].emplace_back(
                new MotorWalkPTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
            else {
                _filRxnTemplates[filType].emplace_back(
                new MotorWalkMTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
            }
        }
        
        //set up reaction templates
        for(auto &r: _chemData.agingReactions[filType]) {
            
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
                auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesFilament[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
                auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesFilament[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesFilament[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    
                    //get position of iterator
                    position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
            _filRxnTemplates[filType].emplace_back(new AgingTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
        }
        
        
        //set up reaction templates
        for(auto &r: _chemData.destructionReactions[filType]) {
            
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
                auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                
                if(it != _chemData.speciesPlusEnd[filType].end()) {
                    //get position of iterator
                    int position = distance(_chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
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
                auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                
                if(it != _chemData.speciesMinusEnd[filType].end()) {
                    //get position of iterator
                    int position = distance(_chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
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
                                      [name](tuple<string, int, double, double, string> element) {
                                          return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    productTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
                }
                
                else if(product.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                      [name](tuple<string, int, double, double, double, string, int> element) {
                                          return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    productTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
                }
                else {
                    cout <<
                    "Third species listed in a destruction reaction must be either bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            
            //add reaction
            _filRxnTemplates[filType].emplace_back(new DestructionTemplate(filType, reactantTemplate, productTemplate, get<2>(r)));
        }
        
        //set up reaction templates
        for(auto &r: _chemData.severingReactions[filType]) {
            
            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;
            
            string reactant = get<0>(r);
            //read strings, and look up type
            
            
            // SPECIES MUST BE FILAMENT
            if(reactant.find("FILAMENT") != string::npos) {
                
                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                int position = 0;
                
                if(it != _chemData.speciesFilament[filType].end()) {
                    //get position of iterator
                    position = distance(_chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
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
            _filRxnTemplates[filType].emplace_back(new SeveringTemplate(filType, reactantTemplate, productTemplate, get<1>(r)));
        }
    }
}

void ChemManager::genFilBindingReactions() {
    
    auto grid = _subSystem->getCompartmentGrid();
    
    //init subsystem ptr
    FilamentBindingManager::_subSystem = _subSystem;
    double rMax, rMin;
    
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        //loop through all compartments
        for(auto C : grid->getCompartments()) {
        
            int managerIndex = 0;
            int linkerIndex = 0;
            int motorIndex = 0;
            
            for(auto &r: _chemData.branchingReactions[filType]) {
                
                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;
                
                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);
                
                cout << reactants.size() << " " << products.size() << endl;
                
                //Checks on number of reactants, products
                if(reactants.size() != BRANCHINGREACTANTS ||
                   products.size() != BRANCHINGPRODUCTS) {
                    cout << "Invalid branching reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                string brancherName;
                
                //FIRST PRODUCT SPECIES MUST BE BRANCHER
                short brancherInt;
                
                auto product = products[0];
                if(product.find("BRANCHER") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(_chemData.speciesBrancher[filType].begin(), _chemData.speciesBrancher[filType].end(), name);
                    
                    if(it != _chemData.speciesBrancher[filType].end()) {
                        
                        brancherName = name;
                        
                        //get position of iterator
                        brancherInt = distance(_chemData.speciesBrancher[filType].begin(), it);
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
                short plusEnd;
                
                product = products[1];
                if(product.find("PLUSEND") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                    
                    if(it != _chemData.speciesPlusEnd[filType].end()) {
                        
                        //get position of iterator
                        plusEnd = distance(_chemData.speciesPlusEnd[filType].begin(), it);
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
                
                int numDiffusingReactant = 0; // Used in determining volume dependence

                //FIRST AND SECOND REACTANTS MUST BE BULK OR DIFFUSING
                auto reactant = reactants[0];
                
                if(reactant.find("BULK") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    
                    auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                      [name](tuple<string, int, double, double, string> element) {
                                          return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }
                
                else if(reactant.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    
                    auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                      [name](tuple<string, int, double, double, double, string, int> element) {
                                          return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
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
                                      [name](tuple<string, int, double, double, string> element) {
                                          return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }
                
                else if(reactant.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                      [name](tuple<string, int, double, double, double, string, int> element) {
                                          return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
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
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                    int position = 0;
                    
                    if(it != _chemData.speciesBound[filType].end()) {
                        
                        //get position of iterator
                        position = distance(_chemData.speciesBound[filType].begin(), it);
                        
                        if(position != SysParams::CParams.brancherBoundIndex[filType]) {
                            cout <<
                            "Third species listed in a branching reaction must be the corresponding brancher empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        //find the species single binding, push
                        string bename = SpeciesNamesDB::genBindingName(brancherName, name);
                        
                        reactantSpecies.push_back(C->findSpeciesByName(bename));
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
                
                //Create reaction
                float onRate = get<2>(r);
                float offRate = get<3>(r);
                auto temp=SysParams::BUBBareRate;
                if(temp.size()>0)
                    temp[brancherInt]=offRate;
                else
                    temp.push_back(offRate);
                SysParams::BUBBareRate=temp;
                //get nucleation zone
                string nzstr = get<4>(r);
                NucleationZoneType nucleationZone;
                if(nzstr == "ALL")
                    nucleationZone = NucleationZoneType::ALL;
                else if(nzstr == "BOUNDARY")
                    nucleationZone = NucleationZoneType::BOUNDARY;
                else if(nzstr == "TOPBOUNDARY")
                    nucleationZone = NucleationZoneType::TOPBOUNDARY;
                else {
                    cout << "Nucleation zone type specified in a branching reaction not valid. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                double nucleationDist = get<5>(r);
                
                ReactionBase* rxn = new Reaction<3,0>(reactantSpecies, onRate, false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::BRANCHING);
                
                C->addInternalReaction(rxn);
                
                //create manager
                BranchingManager* bManager = new BranchingManager(rxn, C, brancherInt, brancherName, filType,
                                                                  nucleationZone, nucleationDist);
                C->addFilamentBindingManager(bManager);
                
                bManager->setMIndex(managerIndex++);
                
                //attach callback
                BranchingCallback bcallback(bManager, plusEnd, onRate, offRate, _subSystem);
                ConnectionBlock rcb(rxn->connect(bcallback,false));
            }
            
            for(auto &r: _chemData.linkerReactions[filType]) {
                
                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;
                
                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);
                
                //Checks on number of reactants, products
                if(reactants.size() != LMBINDINGREACTANTS ||
                   products.size() != LMBINDINGPRODUCTS) {
                    cout << "Invalid linker reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                //FIRST TWO SPECIES IN PRODUCTS MUST BE LINKER
                auto product = products[0];
                
                short linkerInt;
                string linkerName;
                
                if(product.find("LINKER") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    linkerName = name;
                    auto it = find(_chemData.speciesLinker[filType].begin(), _chemData.speciesLinker[filType].end(), name);
                    
                    if(it != _chemData.speciesLinker[filType].end()) {
                        
                        linkerName = name;
                        
                        //get position of iterator
                        linkerInt = distance(_chemData.speciesLinker[filType].begin(), it);
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
                    auto it = find(_chemData.speciesLinker[filType].begin(), _chemData.speciesLinker[filType].end(), name);
                    
                    auto name_prod = products[0].substr(0, products[0].find(":"));
                    
                    if(name != name_prod) {
                        cout <<
                        "Linker species in reactants and products of linker reaction must be same. Exiting." <<
                        endl;
                        exit(EXIT_FAILURE);
                    }
                    if(it == _chemData.speciesLinker[filType].end()) {
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
                
                //FIRST TWO REACTANTS SHOULD BE BOUND
                auto reactant = reactants[0];
                if(reactant.find("BOUND") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                    int position = 0;
                    
                    if(it != _chemData.speciesBound[filType].end()) {
                        
                        //get position of iterator
                        position = distance(_chemData.speciesBound[filType].begin(), it);
                        
                        if(position != SysParams::CParams.linkerBoundIndex[filType]) {
                            cout <<
                            "First species listed in a linker reaction must be the corresponding linker empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        //find the species pair binding, push
                        string lname = SpeciesNamesDB::genBindingName(linkerName, name);
                        
                        reactantSpecies.push_back(C->findSpeciesByName(lname));
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
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                    int position = 0;
                    
                    if(it != _chemData.speciesBound[filType].end()) {
                        
                        //get position of iterator
                        position = distance(_chemData.speciesBound[filType].begin(), it);
                        
                        if(position != SysParams::CParams.linkerBoundIndex[filType]) {
                            cout <<
                            "Second species listed in a linker reaction must be the corresponding linker empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    
                    else if(it == _chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    else if(name != reactants[0]) {
                        cout <<
                        "Both bound species listed in a linker reaction must be the same. Exiting."
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
                
                int numDiffusingReactant = 0; // Used in determining volume dependence
                
                //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
                reactant = reactants[2];
                if(reactant.find("BULK") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                      [name](tuple<string, int, double, double, string> element) {
                                          return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }
                
                else if(reactant.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                      [name](tuple<string, int, double, double, double, string, int> element) {
                                          return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout << "Third species listed in a linker reaction must be bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            
                double onRate = get<2>(r);
                double offRate = get<3>(r);
                //aravind 24, June, 2016.
                auto temp=SysParams::LUBBareRate;
                if(temp.size()>0)
                    temp[linkerInt]=offRate;
                else
                    temp.push_back(offRate);
                SysParams::LUBBareRate=temp;
                rMin = get<4>(r);
                rMax = get<5>(r);
                
                ReactionBase* rxn = new Reaction<2,0>(reactantSpecies, onRate, false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::LINKERBINDING);
                
                C->addInternalReaction(rxn);
                
                //create manager
                LinkerBindingManager* lManager = new LinkerBindingManager(rxn, C, linkerInt, linkerName,
                                                                          filType, rMax, rMin);
                C->addFilamentBindingManager(lManager);
                
                lManager->setNLIndex(linkerIndex++);
                lManager->setMIndex(managerIndex++);
                
                //attach callback
                LinkerBindingCallback lcallback(lManager, onRate, offRate, _subSystem);
                ConnectionBlock rcb(rxn->connect(lcallback,false));
            }
            
            for(auto &r: _chemData.motorReactions[filType]) {
                
                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;
                
                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);
                
                //Checks on number of reactants, products
                if(reactants.size() != LMBINDINGREACTANTS ||
                   products.size() != LMBINDINGPRODUCTS) {
                    cout << "Invalid motor reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                //FIRST TWO SPECIES IN PRODUCTS MUST BE MOTOR
                auto product = products[0];
                
                short motorInt;
                string motorName;
                
                if(product.find("MOTOR") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    motorName = name;
                    auto it = find(_chemData.speciesMotor[filType].begin(), _chemData.speciesMotor[filType].end(), name);
                    
                    if(it != _chemData.speciesMotor[filType].end()) {
                        
                        motorName = name;
                        
                        //get position of iterator
                        motorInt = distance(_chemData.speciesMotor[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A motor species that was included in a reaction was not initialized. Exiting." <<
                        endl;
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
                    auto it = find(_chemData.speciesMotor[filType].begin(), _chemData.speciesMotor[filType].end(), name);
                    
                    auto name_prod = products[0].substr(0, products[0].find(":"));
                    
                    if(name != name_prod) {
                        cout <<
                        "Motor species in reactants and products of motor reaction must be same. Exiting." <<
                        endl;
                        exit(EXIT_FAILURE);
                    }
                    if(it == _chemData.speciesMotor[filType].end()) {
                        cout <<
                        "A motor species that was included in a reaction was not initialized. Exiting." <<
                        endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    cout <<
                    "Fifth species listed in a motor reaction must be motor. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                
                //FIRST TWO REACTANTS SHOULD BE BOUND
                auto reactant = reactants[0];
                if(reactant.find("BOUND") != string::npos) {
                    
                    //look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                    int position = 0;
                    
                    if(it != _chemData.speciesBound[filType].end()) {
                        
                        //get position of iterator
                        position = distance(_chemData.speciesBound[filType].begin(), it);
                        
                        if(position != SysParams::CParams.motorBoundIndex[filType]) {
                            cout <<
                            "First species listed in a motor reaction must be the corresponding motor empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        //find the species pair binding, push
                        string mname = SpeciesNamesDB::genBindingName(motorName, name);
                        
                        reactantSpecies.push_back(C->findSpeciesByName(mname));
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
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
                    int position = 0;
                    
                    if(it != _chemData.speciesBound[filType].end()) {
                        //get position of iterator
                        position = distance(_chemData.speciesBound[filType].begin(), it);
                        
                        if(position != SysParams::CParams.motorBoundIndex[filType]) {
                            cout <<
                            "Second species listed in a motor reaction must be the corresponding motor empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    else if(it == _chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    else if(name != reactants[0]) {
                        cout <<
                        "Both bound species listed in a motor reaction must be the same. Exiting."
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

                int numDiffusingReactant = 0; // Used in determining volume dependence
                
                //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
                reactant = reactants[2];
                if(reactant.find("BULK") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                      [name](tuple<string, int, double, double, string> element) {
                                          return get<0>(element) == name ? true : false; });
                    
                    if(it == _chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }
                
                else if(reactant.find("DIFFUSING") != string::npos) {
                    
                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(_chemData.speciesDiffusing.begin(),_chemData.speciesDiffusing.end(),
                                      [name](tuple<string, int, double, double, double, string, int> element) {
                                          return get<0>(element) == name ? true : false; });
                    if(it == _chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout << "Third species listed in a motor reaction must be bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                
                double onRate = get<2>(r);
                double offRate = get<3>(r);
                //aravind June 24, 2016.
                auto temp=SysParams::MUBBareRate;
                if(temp.size()>0)
                    temp[motorInt]=offRate;
                else
                    temp.push_back(offRate);
                SysParams::MUBBareRate=temp;
                rMin = get<4>(r);
                rMax = get<5>(r);
                
                //multiply by num heads to get rate
                ///CHANGED
                ReactionBase* rxn = new Reaction<2,0>(reactantSpecies, onRate, false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::MOTORBINDING);
                
                C->addInternalReaction(rxn);
                
                //create manager
                MotorBindingManager* mManager = new MotorBindingManager(rxn, C, motorInt, motorName,
                                                                        filType, rMax, rMin);
                C->addFilamentBindingManager(mManager);
                
                mManager->setNLIndex(motorIndex++);
                mManager->setMIndex(managerIndex++);
                
                //attach callback
                MotorBindingCallback mcallback(mManager, onRate, offRate, _subSystem);
                ConnectionBlock rcb(rxn->connect(mcallback,false));
            }
        }
        
        //init neighbor lists
        //get a compartment
        Compartment* C0 = grid->getCompartments()[0];
        
        for(auto &manager : C0->getFilamentBindingManagers()) {
            
            LinkerBindingManager* lManager;
            MotorBindingManager* mManager;
            
            if((lManager = dynamic_cast<LinkerBindingManager*>(manager.get()))) {
                
                auto nl =
                new CylinderCylinderNL(lManager->getRMax() + SysParams::Geometry().cylinderSize[filType],
                                   max(lManager->getRMin() - SysParams::Geometry().cylinderSize[filType], 0.0), true);
                
                //add to subsystem and manager
                LinkerBindingManager::_neighborLists.push_back(nl);
                _subSystem->addNeighborList(nl);
            }
            
            else if((mManager = dynamic_cast<MotorBindingManager*>(manager.get()))) {
                
                auto nl =
                new CylinderCylinderNL(mManager->getRMax() + SysParams::Geometry().cylinderSize[filType],
                                   max(mManager->getRMin() - SysParams::Geometry().cylinderSize[filType], 0.0), true);
                
                //add to subsystem and manager
                MotorBindingManager::_neighborLists.push_back(nl);
                _subSystem->addNeighborList(nl);
            }
        }
    }
}


void ChemManager::genSpecies(Compartment& protoCompartment) {
    
    auto grid = _subSystem->getCompartmentGrid();
    
    // add diffusing species (zero copy number for now)
    for(auto &sd : _chemData.speciesDiffusing) {
        
        auto name = get<0>(sd);
        auto diffRate = get<2>(sd);
        auto rtypeStr = get<5>(sd);
        auto numEvents = get<6>(sd);
        
        RSpeciesType type;
        string rsptype(rtypeStr);
        
        if(rsptype == "AVG")
            type = RSpeciesType::AVG;
        else if (rsptype == "REG")
            type = RSpeciesType::REG;
        
        Species* s = protoCompartment.addSpeciesDiffusing(name, 0, max_ulim, type);
        
        //set num events if averaging
        if(rsptype == "AVG")
            ((RSpeciesAvg*)&s->getRSpecies())->setNumEvents(numEvents);
        
        protoCompartment.setDiffusionRate(name, diffRate);
    }
    
    // add bulk species (zero copy number for now)
    for(auto &sb : _chemData.speciesBulk) {
        
        auto name = get<0>(sb);
        auto rtypeStr = get<4>(sb);
        
        RSpeciesType type;
        string rsptype(rtypeStr);
        
        if(rsptype == "CONST")
            type = RSpeciesType::CONST;
        else if (rsptype == "REG")
            type = RSpeciesType::REG;
        
        grid->addSpeciesBulk(name, 0, max_ulim, type);
    }
    
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        // create single binding and pair binding species
        for(auto &sb : _chemData.speciesBrancher[filType]) {
            
            //look at brancher reaction that is associated with this species
            for(auto &rb : _chemData.branchingReactions[filType]) {
                
                auto reactants = get<0>(rb);
                auto products = get<1>(rb);

                auto sb_bound = products[0].substr(0, products[0].find(":"));
                
                cout << reactants.size() << " " << products.size() << endl;
                
                //basic check because we have not yet checked reactions
                if(reactants.size() != BRANCHINGREACTANTS ||
                   products.size() != BRANCHINGPRODUCTS) {
                    cout << "Invalid branching reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                if(sb_bound == sb) {
                    //look at bound species associated
                    string bound = reactants[2].substr(0, reactants[2].find(":"));
                    
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), bound);
                    
                    //quick check for validity
                    if(it == _chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    //add a single binding species with name sb + bound
                    protoCompartment.addSpeciesSingleBinding(SpeciesNamesDB::genBindingName(sb, bound));
                }
            }
        }
        
        for(auto &sl : _chemData.speciesLinker[filType]) {
            
            //look at linker reaction that is associated with this species
            for(auto &rl : _chemData.linkerReactions[filType]) {
                
                auto reactants = get<0>(rl);
                auto products = get<1>(rl);
                
                //basic check because we have not yet checked reactions
                if(reactants.size() != LMBINDINGREACTANTS ||
                   products.size() != LMBINDINGPRODUCTS) {
                    cout << "Invalid linker reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                auto sl_bound = products[0].substr(0, products[0].find(":"));
                
                if(sl_bound == sl) {
                    
                    //look at bound species associated
                    string bound = reactants[0].substr(0, reactants[0].find(":"));
                    
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), bound);
                    
                    //quick check for validity
                    if(it == _chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    //add a single binding species with name sl + bound
                    protoCompartment.addSpeciesPairBinding(SpeciesNamesDB::genBindingName(sl, bound));
                }
            }
        }
        
        for(auto &sm : _chemData.speciesMotor[filType]) {
            
            //look at brancher reaction that is associated with this species
            for(auto &rm : _chemData.motorReactions[filType]) {
                
                auto reactants = get<0>(rm);
                auto products = get<1>(rm);
                
                //basic check because we have not yet checked reactions
                if(reactants.size() != LMBINDINGREACTANTS ||
                   products.size() != LMBINDINGPRODUCTS) {
                    cout << "Invalid motor reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                
                auto sm_bound = products[0].substr(0, products[0].find(":"));
                
                if(sm_bound == sm) {
                    
                    //look at bound species associated
                    string bound = reactants[0].substr(0, reactants[0].find(":"));
                    
                    auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), bound);
                    
                    //quick check for validity
                    if(it == _chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    //add a single binding species with name sm + bound
                    protoCompartment.addSpeciesPairBinding(SpeciesNamesDB::genBindingName(sm, bound));
                }
            }
        }
    
    }
}

void ChemManager::updateCopyNumbers() {
    
    auto grid = _subSystem->getCompartmentGrid();
    
    //look at copy number for each species
    for(auto &s : _chemData.speciesDiffusing) {
        
        auto name = get<0>(s);
        auto copyNumber = get<1>(s);
        auto releaseTime = get<3>(s);
        auto removalTime = get<4>(s);
        
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
                species->updateReactantPropensities();
                
                copyNumber--;
                
                //set zero copy number
                if(copyNumber == 0) get<1>(s) = 0;
            }
        }
        if(tau() >= removalTime && !areEqual(removalTime,0.0) && get<1>(s) != -1) {
            
            ///remove species from all compartments
            for(auto C : grid->getCompartments()) {
            
                Species* species = C->findSpeciesByName(name);
                
                while(species->getN() > 0) {
                
                    species->down();
                    species->updateReactantPropensities();
                }
            }
            ///set as removed by marking copy number to -1
            get<1>(s) = -1;
        }
    }
    
    for(auto &s : _chemData.speciesBulk) {
        
        auto name = get<0>(s);
        auto copyNumber = get<1>(s);
        auto releaseTime = get<2>(s);
        auto removalTime = get<3>(s);
        
        if(tau() >= releaseTime && copyNumber != 0) {
            
            //find the species, set copy number
            Species* species = grid->findSpeciesBulkByName(name);
            species->setN(copyNumber);
            
            //activate reactions
            species->activateReactantReactions();
            
            //set zero copy number
            get<1>(s) = 0;
        }
        
        if(tau() >= removalTime && !areEqual(removalTime,0.0) && get<1>(s) != -1) {
            
            Species* species = grid->findSpeciesBulkByName(name);
            
            species->setN(0);
            
            //activate reactions
            species->passivateReactantReactions();
            
            ///set as removed by marking copy number to -1
            get<1>(s) = -1;
        }
    }
}

void ChemManager::genGeneralReactions(Compartment& protoCompartment) {
    
    auto grid = _subSystem->getCompartmentGrid();
    
    //go through reactions, add each
    for(auto &r: _chemData.genReactions) {
        
        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;
        
        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);
        
        int numDiffusingReactant = 0; // Used in determining volume dependence
        
        for(auto &reactant : reactants) {
            
            if(reactant.find("BULK") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                  [name](tuple<string, int, double, double, string> element) {
                                   return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
            }
            
            else if(reactant.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it =
                find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                        [name](tuple<string, int, double, double, double, string, int> element) {
                            return get<0>(element) == name ? true : false; });
                if(it == _chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(protoCompartment.findSpeciesByName(name));

                ++numDiffusingReactant;
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
                                  [name](tuple<string, int, double, double, string> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(grid->findSpeciesBulkByName(name));
            }
            
            else if(product.find("DIFFUSING") != string::npos) {
                
                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it =
                find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                        [name](tuple<string, int, double, double, double, string, int> element) {
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
            rxn = new Reaction<1,1>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        else {
            cout <<
            "General reaction specified does not match any existing templates. Exiting."
            <<endl;
            exit(EXIT_FAILURE);
        }
        
        //add to compartment
        protoCompartment.addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void ChemManager::genBulkReactions() {
    
    auto grid = _subSystem->getCompartmentGrid();
    
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
                                  [name](tuple<string, int, double, double, string> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
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
                                  [name](tuple<string, int, double, double, string> element) {
                                      return get<0>(element) == name ? true : false; });
                
                if(it == _chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(grid->findSpeciesBulkByName(name));
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
        grid->addBulkReactionUnique(unique_ptr<ReactionBase>(rxn));
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void ChemManager::genNucleationReactions() {
    
    auto grid = _subSystem->getCompartmentGrid();
    
#if !defined(REACTION_SIGNALING)
    if(!_chemData.nucleationReactions.empty()) {
        
        cout << "Nucleation reactions rely on reaction signaling. Please set this "
        << "compilation macro. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
#endif
    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
    
        //loop through all compartments
        for(auto C : grid->getCompartments()) {
            
            //go through reactions, add each
            for(auto &r: _chemData.nucleationReactions[filType]) {
                
                vector<Species*> reactantSpecies;
                
                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);
                
                if(reactants.size() != NUCLEATIONREACTANTS ||
                   products.size() != NUCLEATIONPRODUCTS) {
                    cout << "Invalid nucleation reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                bool diffusing = false;
                
                int numDiffusingReactant = 0; // Used in determining volume dependence

                for(auto &reactant : reactants) {
                    if(reactant.find("BULK") != string::npos) {
                        
                        //Look up species, make sure in list
                        string name = reactant.substr(0, reactant.find(":"));
                        auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
                                          [name](tuple<string, int, double, double, string> element) {
                                              return get<0>(element) == name ? true : false; });
                        
                        if(it == _chemData.speciesBulk.end()) {
                            cout <<
                            "A bulk species that was included in a reaction was not initialized. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        reactantSpecies.push_back(
                                                  grid->findSpeciesBulkByName(name));
                    }
                    
                    else if(reactant.find("DIFFUSING") != string::npos) {
                        
                        //Look up species, make sure in list
                        string name = reactant.substr(0, reactant.find(":"));
                        auto it =
                        find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
                                [name](tuple<string, int, double, double, double, string, int> element) {
                                    return get<0>(element) == name ? true : false; });
                        if(it == _chemData.speciesDiffusing.end()) {
                            cout <<
                            "A diffusing species that was included in a reaction was not initialized. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        reactantSpecies.push_back(C->findSpeciesByName(name));
                        diffusing = true;

                        ++numDiffusingReactant;
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
                ReactionBase* rxn = new Reaction<2,0>(reactantSpecies, get<2>(r), false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::FILAMENTCREATION);
                
                C->addInternalReaction(rxn);
                
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
                    auto it = find(_chemData.speciesPlusEnd[filType].begin(), _chemData.speciesPlusEnd[filType].end(), name);
                    
                    if(it != _chemData.speciesPlusEnd[filType].end()) {
                        //get position of iterator
                        plusEnd = distance(_chemData.speciesPlusEnd[filType].begin(), it);
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
                    auto it = find(_chemData.speciesFilament[filType].begin(), _chemData.speciesFilament[filType].end(), name);
                    
                    if(it != _chemData.speciesFilament[filType].end()) {
                        //get position of iterator
                        filament = distance(_chemData.speciesFilament[filType].begin(), it);
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
                    auto it = find(_chemData.speciesMinusEnd[filType].begin(), _chemData.speciesMinusEnd[filType].end(), name);
                    
                    if(it != _chemData.speciesMinusEnd[filType].end()) {
                        //get position of iterator
                        minusEnd = distance(_chemData.speciesMinusEnd[filType].begin(), it);
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
                fcallback(plusEnd, filament, minusEnd, filType, _subSystem, creationCompartment);
                ConnectionBlock rcb(rxn->connect(fcallback,false));
    #endif
            }
        }
    }
}

void ChemManager::initializeSystem(ChemSim* chemSim) {
    
    auto grid = _subSystem->getCompartmentGrid();
    
    setupBindingSites();
    configCMonomer();
    
    //Setup all species diffusing and bulk
    Compartment& cProto = grid->getProtoCompartment();
    
    genSpecies(cProto);
    
    //will print reactions as well
    genGeneralReactions(cProto);
    genBulkReactions();
    
    //initialize all compartments equivalent to cproto
    //will copy all general and bulk reactions
    for(auto C : grid->getCompartments())
        *C = cProto;
    for(auto C : grid->getCompartments())
        C->generateAllDiffusionReactions();
    
    //try initial copy number setting
    updateCopyNumbers();
    
    genNucleationReactions();
    genFilBindingReactions();
    
    //add reactions in compartment grid to chemsim
    grid->addChemSimReactions(chemSim);
    
    genFilReactionTemplates();
}


void ChemManager::initializeCCylinder(CCylinder* cc,
                                      bool extensionFront,
                                      bool extensionBack,
                                      bool initialization) {
    
    //get some related objects
    Compartment* C = cc->getCompartment();
    Cylinder* c = cc->getCylinder();
    
    Filament* f = (Filament*)(c->getParent());
    short filType = f->getType();
    
    //add monomers to cylinder
    for(int i = 0; i < cc->getSize(); i++) {
        CMonomer* m = new CMonomer(filType);
        initCMonomer(m, filType, C);
        cc->addCMonomer(m);
        
        if(find(SysParams::Chemistry().bindingSites[filType].begin(),
                SysParams::Chemistry().bindingSites[filType].end(), i)
            !=  SysParams::Chemistry().bindingSites[filType].end()) {
            
            //add callback to all binding sites
            UpdateBrancherBindingCallback bcallback(c, i);
            
            Species* bs = cc->getCMonomer(i)->speciesBound(
                          SysParams::CParams.brancherBoundIndex[filType]);
            ConnectionBlock rcbb(bs->connect(bcallback,false));
            
            UpdateLinkerBindingCallback lcallback(c, i);
            
            Species* ls = cc->getCMonomer(i)->speciesBound(
                          SysParams::CParams.linkerBoundIndex[filType]);
            ConnectionBlock rcbl(ls->connect(lcallback,false));
            
            UpdateMotorBindingCallback mcallback(c, i);
            
            Species* ms = cc->getCMonomer(i)->speciesBound(
                          SysParams::CParams.motorBoundIndex[filType]);
            ConnectionBlock rcbm(ms->connect(mcallback,false));
        }
    }
    
    //get last ccylinder
    CCylinder* lastcc = nullptr;
    
    //extension of front
    if(extensionFront) {
        lastcc = f->getCylinderVector().back()->getCCylinder();
        for(auto &r : _filRxnTemplates[filType]) r->addReaction(lastcc, cc);
    }
    //extension of back
    else if(extensionBack) {
        lastcc = f->getCylinderVector().front()->getCCylinder();
        for(auto &r : _filRxnTemplates[filType]) r->addReaction(cc, lastcc);
    }
    
    //Base case, initialization
    else if (initialization) {
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
            
            for(auto j : SysParams::CParams.bindingIndices[filType])
                m1->speciesBound(j)->up();
            
            //fill new cylinder with default filament value
            for(int i = 0; i < cc->getSize() - 1; i++) {
                cc->getCMonomer(i)->speciesFilament(0)->up();
                
                for(auto j : SysParams::CParams.bindingIndices[filType])
                    cc->getCMonomer(i)->speciesBound(j)->up();
            }
            for(auto &r : _filRxnTemplates[filType]) r->addReaction(lastcc, cc);
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
                
                for(auto j : SysParams::CParams.bindingIndices[filType])
                    cc->getCMonomer(i)->speciesBound(j)->up();
            }
        }
    }
    //Add all reaction templates to this cylinder
    for(auto &r : _filRxnTemplates[filType]) { r->addReaction(cc); }
}

