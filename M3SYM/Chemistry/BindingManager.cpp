
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

BranchingManager::BranchingManager(ReactionBase* reaction, Compartment* compartment,
                                   short boundInt, string boundName)
    : FilamentBindingManager(reaction, compartment, boundInt, boundName) {
    
    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[2]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}

void BranchingManager::updatePossibleBindings(CCylinder* cc, short bindingSite) {
    
    //remove all tuples which have this ccylinder
    int oldN = _bindingSpecies->getN();
    
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));

    //now re add valid binding site
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY)
            _possibleBindings.insert(tuple<CCylinder*, short>(cc, bindingSite));
    
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void BranchingManager::replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) {
    
    vector<tuple<CCylinder*, short>> toAdd;
    
    //remove all tuples which have this ccylinder
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(*it) == oldcc) {
            
            //add new tuple
            toAdd.push_back(tuple<CCylinder*, short>(newcc, get<1>(*it)));
            //erase
            _possibleBindings.erase(it++);
        }
        else ++it;
    }
    //add all to set
    for(auto t : toAdd) _possibleBindings.insert(t);
}

void BranchingManager::updateAllPossibleBindings() {
    
    //clear all
    int oldN = _bindingSpecies->getN();
    
    _possibleBindings.clear();
    
    for(auto &c : _compartment->getCylinders()) {
    
        auto cc = c->getCCylinder();
        
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites.begin();
                 it != SysParams::Chemistry().bindingSites.end(); it++)
            
            if (cc->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY)
                
                _possibleBindings.insert(tuple<CCylinder*, short>(cc, *it));
    }
    
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

LinkerBindingManager::LinkerBindingManager(ReactionBase* reaction, Compartment* compartment,
                                           short boundInt, string boundName, float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[0]->getSpecies().getName();
          
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void LinkerBindingManager::updatePossibleBindings(CCylinder* cc, short bindingSite) {
    
    Cylinder* c = cc->getCylinder();
    
    int oldN = _bindingSpecies->getN();
    
    //remove all tuples which have this ccylinder as key
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    //now re add valid binding sites
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        
        for (auto cn : _nlContainers[_index]->getNeighborList()->getNeighbors(cc->getCylinder())) {
            
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
                    if(c->getID() > cn->getID()) {
                        _possibleBindings.emplace(tuple<CCylinder*, short>(cc, bindingSite),
                                                  tuple<CCylinder*, short>(ccn, *it));
                    }
                    else {
                        if(ccn->getCylinder()->getCompartment() == _compartment) {
                            _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                      tuple<CCylinder*, short>(cc, bindingSite));
                        }
                    }
                }
            }
        }
    }
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void LinkerBindingManager::replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) {
    
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> toAdd;
    
    //remove all tuples which have this ccylinder
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->first) == oldcc) {
            
            //add new tuple
            toAdd.emplace(tuple<CCylinder*, short>(newcc, get<1>(it->first)), it->second);
            //erase
            _possibleBindings.erase(it++);
        }
        else if (get<0>(it->second) == oldcc) {
            
            //add new tuple
            toAdd.emplace(it->first, tuple<CCylinder*, short>(newcc, get<1>(it->second)));
            //erase
            _possibleBindings.erase(it++);
        }
        else ++it;
    }
    
    //add all to set
    for(auto it : toAdd) _possibleBindings.emplace(it.first, it.second);
}

void LinkerBindingManager::updateAllPossibleBindings() {
    
    int oldN = _bindingSpecies->getN();
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
    
        auto cc = c->getCCylinder();
    
        for(auto it1 = SysParams::Chemistry().bindingSites.begin();
                 it1 != SysParams::Chemistry().bindingSites.end(); it1++) {
        
            //now re add valid binding sites
            if (cc->getCMonomer(*it1)->activeSpeciesBound() == BOUND_EMPTY) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _nlContainers[_index]->getNeighborList()->getNeighbors(cc->getCylinder())) {
                    
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
                            else {
                                if(ccn->getCylinder()->getCompartment() == _compartment)
                                    _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it2),
                                                              tuple<CCylinder*, short>(cc, *it1));
                            }
                        }
                    }
                }
            }
        }
    }
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

MotorBindingManager::MotorBindingManager(ReactionBase* reaction, Compartment* compartment,
                                         short boundInt, string boundName, float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[0]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}

void MotorBindingManager::updatePossibleBindings(CCylinder* cc, short bindingSite) {
    
    Cylinder* c = cc->getCylinder();
    
    int oldN = _bindingSpecies->getN();
    
    //remove all tuples which have this ccylinder as key
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    //now re add valid binding sites
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _nlContainers[_index]->getNeighborList()->getNeighbors(cc->getCylinder())) {
            
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
                    if(c->getID() > cn->getID()) {
                        _possibleBindings.emplace(tuple<CCylinder*, short>(cc, bindingSite),
                                                  tuple<CCylinder*, short>(ccn, *it));
                    }
                    else {
                        if(ccn->getCylinder()->getCompartment() == _compartment) {
                            _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it),
                                                      tuple<CCylinder*, short>(cc, bindingSite));
                        }
                    }
                }
            }
        }
    }
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void MotorBindingManager::replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) {
    
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> toAdd;
    
    //remove all tuples which have this ccylinder
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->first) == oldcc) {
            
            //add new tuple
            toAdd.emplace(tuple<CCylinder*, short>(newcc, get<1>(it->first)), it->second);
            //erase
            _possibleBindings.erase(it++);
        }
        else if (get<0>(it->second) == oldcc) {
            
            //add new tuple
            toAdd.emplace(it->first, tuple<CCylinder*, short>(newcc, get<1>(it->second)));
            //erase
            _possibleBindings.erase(it++);
        }
        else ++it;
    }
    
    //add all to set
    for(auto it : toAdd) _possibleBindings.emplace(it.first, it.second);
}


void MotorBindingManager::updateAllPossibleBindings() {
    
    int oldN = _bindingSpecies->getN();
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
        
        auto cc = c->getCCylinder();
        
        for(auto it1 = SysParams::Chemistry().bindingSites.begin();
            it1 != SysParams::Chemistry().bindingSites.end(); it1++) {
            
            //now re add valid binding sites
            if (cc->getCMonomer(*it1)->activeSpeciesBound() == BOUND_EMPTY) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _nlContainers[_index]->getNeighborList()->getNeighbors(cc->getCylinder())) {
                    
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
                            else {
                                if(ccn->getCylinder()->getCompartment() == _compartment)
                                    _possibleBindings.emplace(tuple<CCylinder*, short>(ccn, *it2),
                                                              tuple<CCylinder*, short>(cc, *it1));
                            }
                        }
                    }
                }
            }
        }
    }
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

