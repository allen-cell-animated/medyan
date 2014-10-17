//
//  CMonomer.cpp
//  Cyto
//
//  Created by James Komianos on 9/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CMonomer.h"
#include "SystemParameters.h"

CMonomer::CMonomer() {
    
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    _speciesFilament = new SpeciesFilament*[numFilamentSpecies];
    memset(_speciesFilament, 0, numFilamentSpecies);
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    _speciesPlusEnd = new SpeciesPlusEnd*[numPlusEndSpecies];
    memset(_speciesPlusEnd, 0, numPlusEndSpecies);
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    _speciesMinusEnd = new SpeciesMinusEnd*[numMinusEndSpecies];
    memset(_speciesMinusEnd, 0, numMinusEndSpecies);
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    _speciesBound = new SpeciesBound*[numBoundSpecies];
    memset(_speciesBound, 0, numBoundSpecies);
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    _speciesLinker = new SpeciesLinker*[numLinkerSpecies];
    memset(_speciesLinker, 0, numLinkerSpecies);
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    _speciesMotor = new SpeciesMotor*[numMotorSpecies];
    memset(_speciesMotor, 0, numMotorSpecies);

};

CMonomer::~CMonomer() noexcept{
    
    delete[] _speciesFilament;
    delete[] _speciesPlusEnd;
    delete[] _speciesMinusEnd;
    delete[] _speciesBound;
    delete[] _speciesLinker;
    delete[] _speciesMotor;
}


CMonomer::CMonomer(const CMonomer& rhs, Compartment* c) {
    
    
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = rhs._speciesFilament[i];
        SpeciesFilament* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesFilament(sNew);
    }
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = rhs._speciesPlusEnd[i];
        SpeciesPlusEnd* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesPlusEnd(sNew);
    }
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = rhs._speciesMinusEnd[i];
        SpeciesMinusEnd* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesMinusEnd(sNew);
    }
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = rhs._speciesBound[i];
        SpeciesBound* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesBound(sNew);
    }
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = rhs._speciesLinker[i];
        SpeciesLinker* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesLinker(sNew);
    }
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = rhs._speciesMotor[i];
        SpeciesMotor* sNew = s->clone();
        c->addSpeciesUnique(std::unique_ptr<Species>(sNew));
        addSpeciesMotor(sNew);
    }
}

///Add a species filament
void CMonomer::addSpeciesFilament(SpeciesFilament* s) {
    
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        if(_speciesFilament[i] == nullptr) {
            _speciesFilament[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add filament species to a monomer.\
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species minus end
void CMonomer::addSpeciesPlusEnd(SpeciesPlusEnd* s) {
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        if(_speciesPlusEnd[i] == nullptr) {
            _speciesPlusEnd[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add plus end species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species plus end
void CMonomer::addSpeciesMinusEnd(SpeciesMinusEnd* s) {
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        if(_speciesMinusEnd[i] == nullptr) {
            _speciesMinusEnd[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add minus end species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species bound
void CMonomer::addSpeciesBound(SpeciesBound* s) {
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        if(_speciesBound[i] == nullptr) {
            _speciesBound[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add bound species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}

///Add a species linker
void CMonomer::addSpeciesLinker(SpeciesLinker* s) {
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        if(_speciesLinker[i] == nullptr) {
            _speciesLinker[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add linker species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}
///Add a species motor
void CMonomer::addSpeciesMotor(SpeciesMotor* s) {
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        if(_speciesMotor[i] == nullptr) {
            _speciesMotor[i] = s;
            return;
        }
    }
    ///return error if we get here
    std::cout << "Could not add motor species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << std::endl;
    exit(EXIT_FAILURE);
}


void CMonomer::print()
{
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = _speciesPlusEnd[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = _speciesMinusEnd[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = _speciesLinker[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = _speciesMotor[i];
        if(s != nullptr && s->getN() >= 1) std::cout << s->getName().at(0);
    }
    
}
