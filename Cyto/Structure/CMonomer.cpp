//
//  CMonomer.cpp
//  Cyto
//
//  Created by James Komianos on 9/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CMonomer.h"
CMonomer::CMonomer() {
    memset(_speciesFilament, 0, sizeof(_speciesFilament));
    memset(_speciesPlusEnd, 0, sizeof(_speciesPlusEnd));
    memset(_speciesMinusEnd, 0, sizeof(_speciesMinusEnd));
    memset(_speciesBound, 0, sizeof(_speciesBound));
    memset(_speciesLinker, 0, sizeof(_speciesLinker));
    memset(_speciesMotor, 0, sizeof(_speciesMotor));
};


CMonomer::CMonomer(const CMonomer& rhs, Compartment* c) {
    
    for(int i = 0; i < NUMSPECIESFILAMENT; i++) {
        SpeciesFilament* s = rhs._speciesFilament[i];
        if(s != nullptr) {
            SpeciesFilament* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesFilament(sNew);
        }
    }
    
    for(int i = 0; i < NUMSPECIESPLUSEND; i++) {
        SpeciesPlusEnd* s = rhs._speciesPlusEnd[i];
        if(s != nullptr) {
            SpeciesPlusEnd* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesPlusEnd(sNew);
        }
    }
    
    for(int i = 0; i < NUMSPECIESMINUSEND; i++) {
        SpeciesMinusEnd* s = rhs._speciesMinusEnd[i];
        if(s != nullptr) {
            SpeciesMinusEnd* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesMinusEnd(sNew);
        }
    }
    
    for(int i = 0; i < NUMSPECIESBOUND; i++) {
        SpeciesBound* s = rhs._speciesBound[i];
        if(s != nullptr) {
            SpeciesBound* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesBound(sNew);
        }
    }
    
    for(int i = 0; i < NUMSPECIESLINKER; i++) {
        SpeciesLinker* s = rhs._speciesLinker[i];
        if(s != nullptr) {
            SpeciesLinker* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesLinker(sNew);
        }
    }
    
    for(int i = 0; i < NUMSPECIESMOTOR; i++) {
        SpeciesMotor* s = rhs._speciesMotor[i];
        if(s != nullptr) {
            SpeciesMotor* sNew = s->clone();
            c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
            addSpeciesMotor(sNew);
        }
    }

}

///Add a species filament
void CMonomer::addSpeciesFilament(SpeciesFilament* s) {
    for(int i = 0; i < NUMSPECIESFILAMENT; i++) {
        if(_speciesFilament[i] == nullptr) {
            _speciesFilament[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add filament species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species minus end
void CMonomer::addSpeciesPlusEnd(SpeciesPlusEnd* s) {
    for(int i = 0; i < NUMSPECIESPLUSEND; i++) {
        if(_speciesPlusEnd[i] == nullptr) {
            _speciesPlusEnd[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add plus end species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species plus end
void CMonomer::addSpeciesMinusEnd(SpeciesMinusEnd* s) {
    for(int i = 0; i < NUMSPECIESMINUSEND; i++) {
        if(_speciesMinusEnd[i] == nullptr) {
            _speciesMinusEnd[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add minus end species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species bound
void CMonomer::addSpeciesBound(SpeciesBound* s) {
    
    for(int i = 0; i < NUMSPECIESBOUND; i++) {
        if(_speciesBound[i] == nullptr) {
            _speciesBound[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add bound species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species linker
void CMonomer::addSpeciesLinker(SpeciesLinker* s) {
    for(int i = 0; i < NUMSPECIESLINKER; i++) {
        if(_speciesLinker[i] == nullptr) {
            _speciesLinker[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add linker species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
}
///Add a species motor
void CMonomer::addSpeciesMotor(SpeciesMotor* s) {
    for(int i = 0; i < NUMSPECIESMOTOR; i++) {
        if(_speciesMotor[i] == nullptr) {
            _speciesMotor[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add motor species to a monomer. Exiting" << std::endl;
    exit(EXIT_FAILURE);
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