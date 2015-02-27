
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

#include "CMonomer.h"

#include "CBound.h"

#include "SysParams.h"

CMonomer::CMonomer() {
    
    short numFilamentSpecies =
        SysParams::Chemistry().numFilamentSpecies;
    _speciesFilament = new SpeciesFilament*[numFilamentSpecies]();
    
    short numPlusEndSpecies =
        SysParams::Chemistry().numPlusEndSpecies;
    _speciesPlusEnd = new SpeciesPlusEnd*[numPlusEndSpecies]();
    
    short numMinusEndSpecies =
        SysParams::Chemistry().numMinusEndSpecies;
    _speciesMinusEnd = new SpeciesMinusEnd*[numMinusEndSpecies]();
    
    short numBoundSpecies =
        SysParams::Chemistry().numBoundSpecies;
    _speciesBound = new SpeciesBound*[numBoundSpecies]();
    
    short numLinkerSpecies =
        SysParams::Chemistry().numLinkerSpecies;
    _speciesLinker = new SpeciesLinker*[numLinkerSpecies]();
    
    short numMotorSpecies =
        SysParams::Chemistry().numMotorSpecies;
    _speciesMotor = new SpeciesMotor*[numMotorSpecies]();
    
    short numBrancherSpecies =
        SysParams::Chemistry().numBrancherSpecies;
    _speciesBrancher = new SpeciesBrancher*[numBrancherSpecies]();

};

CMonomer::~CMonomer() noexcept{
    
    delete[] _speciesFilament;
    delete[] _speciesPlusEnd;
    delete[] _speciesMinusEnd;
    delete[] _speciesBound;
    delete[] _speciesLinker;
    delete[] _speciesMotor;
    delete[] _speciesBrancher;
}


CMonomer::CMonomer(const CMonomer& rhs, Compartment* c) : CMonomer() {

    short numFilamentSpecies = SysParams::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = rhs._speciesFilament[i];
        SpeciesFilament* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesFilament(sNew);
    }
    
    short numPlusEndSpecies = SysParams::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = rhs._speciesPlusEnd[i];
        SpeciesPlusEnd* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesPlusEnd(sNew);
    }
    
    short numMinusEndSpecies = SysParams::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = rhs._speciesMinusEnd[i];
        SpeciesMinusEnd* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesMinusEnd(sNew);
    }
    
    //For bound species, transfer the CBound (if any)
    short numBoundSpecies = SysParams::Chemistry().numBoundSpecies;
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
    
    short numLinkerSpecies = SysParams::Chemistry().numLinkerSpecies;
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
    
    short numMotorSpecies = SysParams::Chemistry().numMotorSpecies;
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
    short numBrancherSpecies = SysParams::Chemistry().numBrancherSpecies;
    for(int i = 0; i < numBrancherSpecies; i++) {
        SpeciesBrancher* s = rhs._speciesBrancher[i];
        SpeciesBrancher* sNew = s->clone();
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        addSpeciesBrancher(sNew);
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            if(cBound->getFirstSpecies() == s) cBound->setFirstSpecies(sNew);
            else cBound->setSecondSpecies(sNew);
        }
    }
}

void CMonomer::addSpeciesFilament(SpeciesFilament* s) {
    
    short numFilamentSpecies = SysParams::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        if(_speciesFilament[i] == 0) {
            _speciesFilament[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesPlusEnd(SpeciesPlusEnd* s) {
    
    short numPlusEndSpecies = SysParams::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        if(_speciesPlusEnd[i] == 0) {
            _speciesPlusEnd[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesMinusEnd(SpeciesMinusEnd* s) {
    
    short numMinusEndSpecies = SysParams::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        if(_speciesMinusEnd[i] == 0) {
            _speciesMinusEnd[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesBound(SpeciesBound* s) {
    
    short numBoundSpecies = SysParams::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        if(_speciesBound[i] == 0) {
            _speciesBound[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesLinker(SpeciesLinker* s) {
    
    short numLinkerSpecies = SysParams::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        if(_speciesLinker[i] == 0) {
            _speciesLinker[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesMotor(SpeciesMotor* s) {
    
    short numMotorSpecies= SysParams::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        if(_speciesMotor[i] == 0) {
            _speciesMotor[i] = s;
            return;
        }
    }
}

void CMonomer::addSpeciesBrancher(SpeciesBrancher* s) {
    
    short numBrancherSpecies= SysParams::Chemistry().numBrancherSpecies;
    for(int i = 0; i < numBrancherSpecies; i++) {
        if(_speciesBrancher[i] == 0) {
            _speciesBrancher[i] = s;
            return;
        }
    }
}


void CMonomer::print()
{
    short numFilamentSpecies = SysParams::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numPlusEndSpecies = SysParams::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = _speciesPlusEnd[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numMinusEndSpecies = SysParams::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = _speciesMinusEnd[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numBoundSpecies = SysParams::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numLinkerSpecies = SysParams::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = _speciesLinker[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numMotorSpecies = SysParams::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = _speciesMotor[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
    short numBrancherSpecies = SysParams::Chemistry().numBrancherSpecies;
    for(int i = 0; i < numBrancherSpecies; i++) {
        SpeciesBrancher* s = _speciesBrancher[i];
        if(s != nullptr && s->getN() >= 1) cout << s->getName();
    }
    
}

short CMonomer::activeSpeciesFilament() {
    short numFilamentSpecies = SysParams::Chemistry().numFilamentSpecies;
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesPlusEnd() {
    short numPlusEndSpecies = SysParams::Chemistry().numPlusEndSpecies;
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesPlusEnd* s = _speciesPlusEnd[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesMinusEnd() {
    short numMinusEndSpecies = SysParams::Chemistry().numMinusEndSpecies;
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesMinusEnd* s = _speciesMinusEnd[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}

short CMonomer::activeSpeciesBound() {
    short numBoundSpecies = SysParams::Chemistry().numBoundSpecies;
    for(int i = 0; i < numBoundSpecies; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesLinker() {
    short numLinkerSpecies = SysParams::Chemistry().numLinkerSpecies;
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesLinker* s = _speciesLinker[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
    
}
short CMonomer::activeSpeciesMotor() {
    short numMotorSpecies = SysParams::Chemistry().numMotorSpecies;
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesMotor* s = _speciesMotor[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesBrancher() {
    short numBrancherSpecies = SysParams::Chemistry().numBrancherSpecies;
    for(int i = 0; i < numBrancherSpecies; i++) {
        SpeciesBrancher* s = _speciesBrancher[i];
        if(s != nullptr && s->getN() >= 1) return i;
    }
    return -1;
}



