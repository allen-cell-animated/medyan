
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

#include "CMonomer.h"

#include "CBound.h"

#include "SystemParameters.h"

CMonomer::CMonomer() {
    
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    _speciesFilament = new SpeciesFilament*[numFilamentSpecies]();
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    _speciesPlusEnd = new SpeciesPlusEnd*[numPlusEndSpecies]();
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    _speciesMinusEnd = new SpeciesMinusEnd*[numMinusEndSpecies]();
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    _speciesBound = new SpeciesBound*[numBoundSpecies]();
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    _speciesLinker = new SpeciesLinker*[numLinkerSpecies]();
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    _speciesMotor = new SpeciesMotor*[numMotorSpecies]();

};

CMonomer::~CMonomer() noexcept{
    
    delete[] _speciesFilament;
    delete[] _speciesPlusEnd;
    delete[] _speciesMinusEnd;
    delete[] _speciesBound;
    delete[] _speciesLinker;
    delete[] _speciesMotor;
}


CMonomer::CMonomer(const CMonomer& rhs, Compartment* c) : CMonomer() {

    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = rhs._speciesFilament[i];
        SpeciesFilament* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesFilament(sNew);
    }
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = rhs._speciesPlusEnd[i];
        SpeciesPlusEnd* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesPlusEnd(sNew);
    }
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = rhs._speciesMinusEnd[i];
        SpeciesMinusEnd* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesMinusEnd(sNew);
    }
    
    //For bound species, transfer the CBound (if any)
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = rhs._speciesBound[i];
        SpeciesBound* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesBound(sNew);
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            if(cBound->getFirstSpecies() == s) cBound->setFirstSpecies(sNew);
            else cBound->setSecondSpecies(sNew);
        }
        
    }
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = rhs._speciesLinker[i];
        SpeciesLinker* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesLinker(sNew);
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            if(cBound->getFirstSpecies() == s) cBound->setFirstSpecies(sNew);
            else cBound->setSecondSpecies(sNew);
        }
    }
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = rhs._speciesMotor[i];
        SpeciesMotor* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesMotor(sNew);
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            if(cBound->getFirstSpecies() == s) cBound->setFirstSpecies(sNew);
            else cBound->setSecondSpecies(sNew);
        }
    }
}

void CMonomer::addSpeciesFilament(SpeciesFilament* s) {
    
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        if(_speciesFilament[i] == 0) {
            _speciesFilament[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add filament species to a monomer.\
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}

void CMonomer::addSpeciesPlusEnd(SpeciesPlusEnd* s) {
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        if(_speciesPlusEnd[i] == 0) {
            _speciesPlusEnd[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add plus end species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}

void CMonomer::addSpeciesMinusEnd(SpeciesMinusEnd* s) {
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        if(_speciesMinusEnd[i] == 0) {
            _speciesMinusEnd[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add minus end species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}

void CMonomer::addSpeciesBound(SpeciesBound* s) {
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        if(_speciesBound[i] == 0) {
            _speciesBound[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add bound species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}

void CMonomer::addSpeciesLinker(SpeciesLinker* s) {
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        if(_speciesLinker[i] == 0) {
            _speciesLinker[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add linker species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}

void CMonomer::addSpeciesMotor(SpeciesMotor* s) {
    
    short numMotorSpecies= SystemParameters::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        if(_speciesMotor[i] == 0) {
            _speciesMotor[i] = s;
            return;
        }
    }
    //return error if we get here
    cout << "Could not add motor species to a monomer. \
        Check that the numer of species in the system input file matches the chemistry input Exiting" << endl;
    exit(EXIT_FAILURE);
}


void CMonomer::print()
{
    short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
    short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = _speciesPlusEnd[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
    short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = _speciesMinusEnd[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
    short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
    short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = _speciesLinker[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
    short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = _speciesMotor[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName().at(0);
    }
    
}
