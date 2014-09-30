//
//  CMonomer.cpp
//  Cyto
//
//  Created by James Komianos on 9/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CMonomer.h"

CMonomer::CMonomer(const CMonomer& rhs, Compartment* c) {
    
    for(auto &s: rhs._speciesFilament) {
        SpeciesFilament* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        _speciesFilament.push_back(sNew);
    }
    for(auto &s: rhs._speciesBound) {
        SpeciesBound* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        _speciesBound.push_back(sNew);
    }
    for(auto &s: rhs._speciesPlusEnd) {
        SpeciesPlusEnd* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        _speciesPlusEnd.push_back(sNew);
    }
    for(auto &s: rhs._speciesMinusEnd) {
        SpeciesMinusEnd* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        _speciesMinusEnd.push_back(sNew);
    }
}

void CMonomer::print()
{
    for (auto &s : _speciesFilament)
        if(s->getN() >= 1) std::cout << s->getName().at(0);
    for (auto &s : _speciesBound)
        if(s->getN() >= 1) std::cout << s->getName().at(0);
    for (auto &s : _speciesPlusEnd)
        if(s->getN() >= 1) std::cout << s->getName().at(0);
    for (auto &s : _speciesMinusEnd)
        if(s->getN() >= 1) std::cout << s->getName().at(0);
}