
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

#include "CCylinder.h"

#include "CBound.h"

CCylinder::CCylinder(const CCylinder& rhs, Compartment* c) : _compartment(c), _pCylinder(rhs._pCylinder)
{
    //copy all monomers, bounds
    for(auto &m : rhs._monomers)
        _monomers.push_back(unique_ptr<CMonomer>(m->clone(c)));
    
    //copy all internal reactions
    for(auto &r: rhs._internalReactions)
        addInternalReaction(r->clone(c->getSpeciesContainer()));
    
    //copy all cross-cylinder reactions
    for(auto it = rhs._crossCylinderReactions.begin();
             it != rhs._crossCylinderReactions.end(); it++) {
        
        for(auto &r : it->second) {
            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            if(r->getCBound() != nullptr) {
                r->getCBound()->setOffReaction(rxnClone);
                rxnClone->setCBound(r->getCBound());
            }
            addCrossCylinderReaction(it->first, rxnClone);
        }
    }
    
    //Copy reacting cylinders, Clone reactions where this cylinder is involved
    for(auto &ccyl : rhs._reactingCylinders) {

        //clone reactions
        for(auto &r: ccyl->getCrossCylinderReactions()[const_cast<CCylinder*>(&rhs)]) {
            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            if(r->getCBound() != nullptr) {
                r->getCBound()->setOffReaction(rxnClone);
                rxnClone->setCBound(r->getCBound());
            }
            ccyl->addCrossCylinderReaction(this, rxnClone);
        }
    }
}

void CCylinder::addInternalReaction(ReactionBase* r) {
    //add to compartment and chemsim
    _compartment->addInternalReactionUnique(unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(r);
    
    //add to local reaction list
    _internalReactions.insert(r);
    
    //activate reaction
    r->activateReaction();
}

void CCylinder::addCrossCylinderReaction(CCylinder* other, ReactionBase* r) {
    //add to compartment and chemsim
    _compartment->addInternalReactionUnique(unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(r);
    
    //add to this reaction map
    _crossCylinderReactions[other].insert(r);
    other->addReactingCylinder(this);
    
    //activate reaction
    r->activateReaction();
}

void CCylinder::addReactingCylinder(CCylinder* other) {
    _reactingCylinders.insert(other);
}


void CCylinder::removeInternalReaction(ReactionBase* r) {
    
    //remove from internal reaction list
    if (_internalReactions.find(r) != _internalReactions.end()) {
        
        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(r);
            
        _internalReactions.erase(r);
    }
}

void CCylinder:: removeAllInternalReactions() {
    
    auto tempReactions = _internalReactions;
    for (auto &r : tempReactions) removeInternalReaction(r);
}

void CCylinder::removeCrossCylinderReaction(CCylinder* other, ReactionBase* r) {
    
    if(_crossCylinderReactions[other].find(r) != _crossCylinderReactions[other].end()) {
       
        _crossCylinderReactions[other].erase(r);
        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(r);
        
        //if number of reactions in cross-cylinder has dropped to zero, delete it
        if(_crossCylinderReactions[other].size() == 0) {
            
            _crossCylinderReactions.erase(other);
            
            //also remove from reacting list of other ccylinder
            auto it = other->_reactingCylinders.find(this);
            if(it != other->_reactingCylinders.end())
                other->_reactingCylinders.erase(it);
        }
    }
}

void CCylinder::removeCrossCylinderReactions(CCylinder* other, bool bindingOnly) {
    
    auto tempReactions = _crossCylinderReactions[other];
    for(auto &r : tempReactions) {
        
        if(bindingOnly) {
            if (!(r->getReactionType() == ReactionType::LINKERUNBINDING
                || r->getReactionType() == ReactionType::MOTORUNBINDING))

                removeCrossCylinderReaction(other, r);
        }
        else removeCrossCylinderReaction(other, r);
    }
}

void CCylinder::removeAllCrossCylinderReactions() {
    
    auto tempMap = _crossCylinderReactions;
    for(auto it = tempMap.begin(); it != tempMap.end(); it++)

        removeCrossCylinderReactions(it->first);
}

void CCylinder::removeReactingCylinder(CCylinder* other) {
    other->removeCrossCylinderReactions(this);
}

void CCylinder::removeAllReactingCylinders() {
    auto tempReactingCylinders = _reactingCylinders;
    for(auto &cc : tempReactingCylinders)
        cc->removeCrossCylinderReactions(this);
}

CCylinder::~CCylinder() {
    //Remove all reactions owned by this ccylinder
    removeAllInternalReactions();
    removeAllCrossCylinderReactions();
    
    //remove all reactions involving this ccylinder
    removeAllReactingCylinders();
    
    //Remove all species
    for(auto &m: _monomers) {
        
        short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
        for(int i = 0; i < numFilamentSpecies; i++) {
            SpeciesFilament* s = m->speciesFilament(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
        short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
        for(int i = 0; i < numPlusEndSpecies; i++) {
            SpeciesPlusEnd* s = m->speciesPlusEnd(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
        short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
        for(int i = 0; i < numMinusEndSpecies; i++) {
            SpeciesMinusEnd* s = m->speciesMinusEnd(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
        short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
        for(int i = 0; i < numBoundSpecies; i++) {
            SpeciesBound* s = m->speciesBound(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
        short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
        for(int i = 0; i < numLinkerSpecies; i++) {
            SpeciesLinker* s = m->speciesLinker(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
        short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
        for(int i = 0; i < numMotorSpecies; i++) {
            SpeciesMotor* s = m->speciesMotor(i);
            if(s != nullptr) {_compartment->removeSpecies(s);}
        }
    }
}

void CCylinder::printCCylinder()
{
    cout << "Compartment:" << _compartment << endl;
    
    cout << "Composition of CCylinder: " << endl;
    for (auto &m : _monomers){
        m->print();
        cout << ":";
    }
    cout << endl;
}


