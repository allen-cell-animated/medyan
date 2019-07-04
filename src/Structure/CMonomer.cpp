//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CMonomer.h"

#include "CBound.h"

#include "SysParams.h"

CMonomer::CMonomer(short filamentType) : _filamentType(filamentType) {
    
    _speciesFilament = new SpeciesFilament*[_numFSpecies[_filamentType]]();
    _speciesBound = new SpeciesBound*[_numBSpecies[_filamentType]]();
};

CMonomer::~CMonomer() noexcept{
    
    delete[] _speciesFilament;
    delete[] _speciesBound;
}


CMonomer::CMonomer(const CMonomer& rhs, Compartment* c)

    : CMonomer(rhs._filamentType) {

    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        
        //clone and add to array
        SpeciesFilament* s = rhs._speciesFilament[i];
        SpeciesFilament* sNew = s->clone();
        
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        _speciesFilament[i] = sNew;
    }
    
    //For bound species, transfer the CBound (if any)
    
    for(int i = 0; i < _numBSpecies[_filamentType]; i++) {
        
        //clone and add to array
        SpeciesBound* s = rhs._speciesBound[i];
        SpeciesBound* sNew = s->clone();
        
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        _speciesBound[i] = sNew;
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            //set species
            if(cBound->getFirstSpecies() == s)
                cBound->setFirstSpecies(sNew);
            else
                cBound->setSecondSpecies(sNew);
        }
    }
}

//PRINT

void CMonomer::print()
{
    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        SpeciesFilament* s = _speciesFilament[i];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            cout << s->getName();
    }
    for(int i = 0; i < _numBSpecies[_filamentType]; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            cout << s->getName();
    }
}

//GETTERS

SpeciesFilament* CMonomer::speciesFilament(int index) {
    short offset = _speciesFilamentIndex[_filamentType][SPECIESFILAMENT];
    return _speciesFilament[index + offset];
}
SpeciesFilament* CMonomer::speciesPlusEnd (int index) {
    short offset = _speciesFilamentIndex[_filamentType][SPECIESPLUSEND];
    return _speciesFilament[index + offset];
}
SpeciesFilament* CMonomer::speciesMinusEnd(int index) {
    short offset = _speciesFilamentIndex[_filamentType][SPECIESMINUSEND];
    return _speciesFilament[index + offset];
}

SpeciesBound* CMonomer::speciesBound(int index) {
    short offset = _speciesBoundIndex[_filamentType][SPECIESBOUND];
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesLinker(int index) {
    short offset = _speciesBoundIndex[_filamentType][SPECIESLINKER];
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesMotor(int index) {
    short offset = _speciesBoundIndex[_filamentType][SPECIESMOTOR];
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesBrancher(int index) {
    short offset = _speciesBoundIndex[_filamentType][SPECIESBRANCHER];
    return _speciesBound[index + offset];
}


//GET ACTIVE

short CMonomer::activeSpeciesFilament() {
    short numFilamentSpecies = SysParams::Chemistry().numFilamentSpecies[_filamentType];
    short offset = _speciesFilamentIndex[_filamentType][SPECIESFILAMENT];
    
    for(int i = 0; i < numFilamentSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesPlusEnd() {
    short numPlusEndSpecies = SysParams::Chemistry().numPlusEndSpecies[_filamentType];
    short offset = _speciesFilamentIndex[_filamentType][SPECIESPLUSEND];
    
    for(int i = 0; i < numPlusEndSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            return i;
    }
    return -1;
}
short CMonomer::activeSpeciesMinusEnd() {
    short numMinusEndSpecies = SysParams::Chemistry().numMinusEndSpecies[_filamentType];
    short offset = _speciesFilamentIndex[_filamentType][SPECIESMINUSEND];
    
    for(int i = 0; i < numMinusEndSpecies; i++) {
        SpeciesFilament* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            return i;
    }
    return -1;
}

short CMonomer::activeSpeciesLinker() {
    short numLinkerSpecies = SysParams::Chemistry().numLinkerSpecies[_filamentType];
    short offset = _speciesBoundIndex[_filamentType][SPECIESLINKER];
    
    for(int i = 0; i < numLinkerSpecies; i++) {
        SpeciesBound* s = _speciesBound[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
    
}
short CMonomer::activeSpeciesMotor() {
    short numMotorSpecies = SysParams::Chemistry().numMotorSpecies[_filamentType];
    short offset = _speciesBoundIndex[_filamentType][SPECIESMOTOR];
    
    for(int i = 0; i < numMotorSpecies; i++) {
        SpeciesBound* s = _speciesBound[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesBrancher() {
    short numBrancherSpecies = SysParams::Chemistry().numBrancherSpecies[_filamentType];
    short offset = _speciesBoundIndex[_filamentType][SPECIESBRANCHER];
    
    for(int i = 0; i < numBrancherSpecies; i++) {
        SpeciesBound* s = _speciesBound[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
}

bool CMonomer::isConsistent() {
    
    //check all species between 0 and 1 inclusive
    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        
        if(!areEqual(_speciesFilament[i]->getN(), 1.0) &&
           !areEqual(_speciesFilament[i]->getN(), 0.0)) {
            
            cout << _speciesFilament[i]->getName() << " has an invalid copy number. It is = "
                 << _speciesFilament[i]->getN() << " and is at species index " << i << "." << endl;
            
            return false;
        }
    }
    //check filament species
    if(activeSpeciesFilament() != -1 &&
       (activeSpeciesPlusEnd() != -1 ||
        activeSpeciesMinusEnd() != -1)) {
           
        cout << "Has a simultaneously active filament and plus/minus end species." << endl;
           
        return false;
    }
    
    return true;
}

vector<vector<short>> CMonomer::_speciesFilamentIndex = vector<vector<short>>(MAX_FILAMENT_TYPES);
vector<vector<short>> CMonomer::_speciesBoundIndex    = vector<vector<short>>(MAX_FILAMENT_TYPES);

vector<short> CMonomer::_numFSpecies = vector<short>(MAX_FILAMENT_TYPES);
vector<short> CMonomer::_numBSpecies = vector<short>(MAX_FILAMENT_TYPES);



