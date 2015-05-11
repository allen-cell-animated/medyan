
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

#include "BindingManager.h"

#include "Compartment.h"
#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

mt19937* FilamentBindingManager::_eng = 0;

vector<CCNLContainer*> LinkerBindingManager::_nlContainers;
vector<CCNLContainer*> MotorBindingManager::_nlContainers;


//BRANCHER

BranchingManager::BranchingManager(ReactionBase* reaction, Compartment* compartment,
                                   short boundInt, string boundName)
    : FilamentBindingManager(reaction, compartment, boundInt, boundName) {
    
    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[2]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}

void BranchingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    //add valid site
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY)
        _possibleBindings.insert(tuple<CCylinder*, short>(cc, bindingSite));
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void BranchingManager::addPossibleBindings(CCylinder* cc) {
    
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
        
        addPossibleBindings(cc, *bit);
}

void BranchingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    //remove tuple which has this ccylinder
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}


void BranchingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
        
        removePossibleBindings(cc, *bit);
}


void BranchingManager::updateAllPossibleBindings() {
    
    //clear all
    _possibleBindings.clear();
    
    for(auto &c : _compartment->getCylinders()) {
    
        auto cc = c->getCCylinder();
        
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites.begin();
                 it != SysParams::Chemistry().bindingSites.end(); it++)
            
            if (cc->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY)
                
                _possibleBindings.insert(tuple<CCylinder*, short>(cc, *it));
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

//LINKER

LinkerBindingManager::LinkerBindingManager(ReactionBase* reaction, Compartment* compartment,
                                           short boundInt, string boundName, float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[0]->getSpecies().getName();
          
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void LinkerBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                       getNeighbors(cc->getCylinder())) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getFilament() == c->getFilament()) continue;
            
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites.begin();
                     it != SysParams::Chemistry().bindingSites.end(); it++) {
                
                if (ccn->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY) {
                    
                    //check distances..
                    auto midpoint1 = (float)bindingSite / SysParams::Geometry().cylinderIntSize;
                    auto pos1 = midPointCoordinate(c->getFirstBead()->coordinate,
                                                   c->getSecondBead()->coordinate, midpoint1);
                    
                    auto midpoint2 = (float)*it / SysParams::Geometry().cylinderIntSize;
                    auto pos2 = midPointCoordinate(cn->getFirstBead()->coordinate,
                                                   cn->getSecondBead()->coordinate, midpoint2);
                    
                    double dist = twoPointDistance(pos1, pos2);
                    
                    if(dist > _rMax || dist < _rMin) continue;
                    
                    //add in correct order
                    if(c->getID() > cn->getID())
                        _possibleBindings.emplace(tuple<CCylinder*, short>(cc, bindingSite),
                                                  tuple<CCylinder*, short>(ccn, *it));
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {
                            
                            _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                      tuple<CCylinder*, short>(cc, bindingSite));
                        }
                        //add in other
                        else {
                            
                            auto m = (LinkerBindingManager*)cn->getCompartment()->
                                      getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);
                            
                            m->_possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                         tuple<CCylinder*, short>(cc, bindingSite));
                        }
                    }
                }
            }
        }
    }
    
    //update affected
    for(auto m : affectedManagers) {
        
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        
        m->updateBindingReaction(oldNOther, newNOther);
        
    }
    
    //update this manager
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void LinkerBindingManager::addPossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
    
        addPossibleBindings(cc, *bit);
}



void LinkerBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;
    
    //remove all tuples which have this ccylinder as key
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    //remove all neighbors which have this binding site pair
    for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                   getNeighbors(cc->getCylinder())) {
        
        if(cn->getCompartment() != _compartment) {
            
            auto m = (LinkerBindingManager*)cn->getCompartment()->
                      getFilamentBindingManagers()[_mIndex].get();
            
            affectedManagers.push_back(m);
        }
    }
    
    //remove, update affected
    for(auto m : affectedManagers) {
        
        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {
            
            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindings.erase(it++);
            
            else ++it;
        }
        
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        
        m->updateBindingReaction(oldNOther, newNOther);
        
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void LinkerBindingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void LinkerBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
    
        auto cc = c->getCCylinder();
    
        for(auto it1 = SysParams::Chemistry().bindingSites.begin();
                 it1 != SysParams::Chemistry().bindingSites.end(); it1++) {
        
            //now re add valid binding sites
            if (cc->getCMonomer(*it1)->activeSpeciesBound() == BOUND_EMPTY) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                               getNeighbors(cc->getCylinder())) {
                    
                    if(cn->getFilament() == c->getFilament()) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites.begin();
                             it2 != SysParams::Chemistry().bindingSites.end(); it2++) {
                        
                        if (ccn->getCMonomer(*it2)->activeSpeciesBound() == BOUND_EMPTY) {
                            
                            //check distances..
                            auto midpoint1 = (float)*it1 / SysParams::Geometry().cylinderIntSize;
                            auto pos1 = midPointCoordinate(c->getFirstBead()->coordinate,
                                                           c->getSecondBead()->coordinate, midpoint1);
                            
                            auto midpoint2 = (float)*it2 / SysParams::Geometry().cylinderIntSize;
                            auto pos2 = midPointCoordinate(cn->getFirstBead()->coordinate,
                                                           cn->getSecondBead()->coordinate, midpoint2);
                            
                            double dist = twoPointDistance(pos1, pos2);
                            
                            if(dist > _rMax || dist < _rMin) continue;
                            
                            //add in correct order
                            if(c->getID() > cn->getID())
                                _possibleBindings.emplace(tuple<CCylinder*, short>(cc, *it1),
                                                          tuple<CCylinder*, short>(ccn, *it2));
                        }
                    }
                }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}


//MOTOR

MotorBindingManager::MotorBindingManager(ReactionBase* reaction, Compartment* compartment,
                                         short boundInt, string boundName, float rMax, float rMin)

: FilamentBindingManager(reaction, compartment, boundInt, boundName),
_rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[0]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void MotorBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                       getNeighbors(cc->getCylinder())) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getFilament() == c->getFilament()) continue;
        
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites.begin();
                     it != SysParams::Chemistry().bindingSites.end(); it++) {
                
                if (ccn->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY) {
                    
                    //check distances..
                    auto midpoint1 = (float)bindingSite / SysParams::Geometry().cylinderIntSize;
                    auto pos1 = midPointCoordinate(c->getFirstBead()->coordinate,
                                                   c->getSecondBead()->coordinate, midpoint1);
                    
                    auto midpoint2 = (float)*it / SysParams::Geometry().cylinderIntSize;
                    auto pos2 = midPointCoordinate(cn->getFirstBead()->coordinate,
                                                   cn->getSecondBead()->coordinate, midpoint2);
                    
                    double dist = twoPointDistance(pos1, pos2);
                    
                    if(dist > _rMax || dist < _rMin) continue;
                    
                    //add in correct order
                    if(c->getID() > cn->getID())
                        _possibleBindings.emplace(tuple<CCylinder*, short>(cc, bindingSite),
                                                  tuple<CCylinder*, short>(ccn, *it));
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {
                            
                            _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                      tuple<CCylinder*, short>(cc, bindingSite));
                        }
                        //add in other
                        else {
                            
                            auto m = (MotorBindingManager*)cn->getCompartment()->
                            getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);
                            
                            m->_possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                         tuple<CCylinder*, short>(cc, bindingSite));
                        }
                    }
                }
            }
        }
    }
    
    //update affected
    for(auto m : affectedManagers) {
        
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        
        m->updateBindingReaction(oldNOther, newNOther);
        
    }
    
    //update this manager
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void MotorBindingManager::addPossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
        
        addPossibleBindings(cc, *bit);
}


void MotorBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
    
    //remove all tuples which have this ccylinder as key
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    //remove all neighbors which have this binding site pair
    for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                   getNeighbors(cc->getCylinder())) {
        
        if(cn->getCompartment() != _compartment) {
            
            auto m = (MotorBindingManager*)cn->getCompartment()->
                      getFilamentBindingManagers()[_mIndex].get();
            
            affectedManagers.push_back(m);
        }
    }
    
    //remove, update affected
    for(auto m : affectedManagers) {
        
        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {
            
            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindings.erase(it++);
            
            else ++it;
        }
        
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        
        m->updateBindingReaction(oldNOther, newNOther);
        
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void MotorBindingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites.begin();
             bit != SysParams::Chemistry().bindingSites.end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void MotorBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
        
        auto cc = c->getCCylinder();
        
        for(auto it1 = SysParams::Chemistry().bindingSites.begin();
                 it1 != SysParams::Chemistry().bindingSites.end(); it1++) {
            
            //now re add valid binding sites
            if (cc->getCMonomer(*it1)->activeSpeciesBound() == BOUND_EMPTY) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _nlContainers[_nlIndex]->getNeighborList()->
                               getNeighbors(cc->getCylinder())) {
                    
                    if(cn->getFilament() == c->getFilament()) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites.begin();
                             it2 != SysParams::Chemistry().bindingSites.end(); it2++) {
                        
                        if (ccn->getCMonomer(*it2)->activeSpeciesBound() == BOUND_EMPTY) {
                            
                            //check distances..
                            auto midpoint1 = (float)*it1 / SysParams::Geometry().cylinderIntSize;
                            auto pos1 = midPointCoordinate(c->getFirstBead()->coordinate,
                                                           c->getSecondBead()->coordinate, midpoint1);
                            
                            auto midpoint2 = (float)*it2 / SysParams::Geometry().cylinderIntSize;
                            auto pos2 = midPointCoordinate(cn->getFirstBead()->coordinate,
                                                           cn->getSecondBead()->coordinate, midpoint2);
                            
                            double dist = twoPointDistance(pos1, pos2);
                            
                            if(dist > _rMax || dist < _rMin) continue;
                            
                            //add in correct order
                            if(c->getID() > cn->getID())
                                _possibleBindings.emplace(tuple<CCylinder*, short>(cc, *it1),
                                                          tuple<CCylinder*, short>(ccn, *it2));
                        }
                    }
                }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

