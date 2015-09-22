
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

#include "SubSystem.h"
#include "Boundary.h"

#include "MathFunctions.h"
#include "GController.h"
#include "SysParams.h"

using namespace mathfunc;

mt19937* FilamentBindingManager::_eng = 0;
SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CCNeighborList*> LinkerBindingManager::_neighborLists;
vector<CCNeighborList*> MotorBindingManager::_neighborLists;

//BRANCHER

BranchingManager::BranchingManager(ReactionBase* reaction,
                                   Compartment* compartment,
                                   short boundInt, string boundName,
                                   NucleationZoneType zone,
                                   double nucleationDistance)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _nucleationZone(zone), _nucleationDistance(nucleationDistance) {
    
    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[SPECIESBB_BINDING_INDEX]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}

void BranchingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    
    bool inZone = true;
    //see if in nucleation zone
    if(_nucleationZone != NucleationZoneType::ALL) {
        
        auto mp = (float)bindingSite / SysParams::Geometry().cylinderIntSize;
        
        auto x1 = cc->getCylinder()->getFirstBead()->coordinate;
        auto x2 = cc->getCylinder()->getSecondBead()->coordinate;
        
        auto coord = midPointCoordinate(x1, x2, mp);
        
        //set nucleation zone
        if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {
            
            //if top boundary, check if we are above the center coordinate in z
            if(_nucleationZone == NucleationZoneType::TOPBOUNDARY) {
                
                if(coord[2] >= GController::getCenter()[2])
                    inZone = true;
                else
                    inZone = false;
            }
            else inZone = true;
        }
        else
            inZone = false;
    }
    
    //add valid site
    if (cc->getCMonomer(bindingSite)->speciesBound(B_BINDING_INDEX)->getN() == 1 && inZone) {
        
        auto t = tuple<CCylinder*, short>(cc, bindingSite);
        _possibleBindings.insert(t);
    }
    
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
                 it != SysParams::Chemistry().bindingSites.end(); it++) {
            
            bool inZone = true;
            //see if in nucleation zone
            if(_nucleationZone != NucleationZoneType::ALL) {
                
                auto mp = (float)*it / SysParams::Geometry().cylinderIntSize;
                
                auto x1 = cc->getCylinder()->getFirstBead()->coordinate;
                auto x2 = cc->getCylinder()->getSecondBead()->coordinate;
                
                auto coord = midPointCoordinate(x1, x2, mp);
                
                //set nucleation zone
                if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {
                    
                    //if top boundary, check if we are above the center coordinate in z
                    if(_nucleationZone == NucleationZoneType::TOPBOUNDARY) {
                        
                        if(coord[2] >= GController::getCenter()[2])
                            inZone = true;
                        else
                            inZone = false;
                    }
                    else inZone = true;
                }
                else
                    inZone = false;
            }
            if (cc->getCMonomer(*it)->speciesBound(B_BINDING_INDEX)->getN() == 1 && inZone) {
                
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindings.insert(t);
            }
        }
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

//LINKER

LinkerBindingManager::LinkerBindingManager(ReactionBase* reaction,
                                           Compartment* compartment,
                                           short boundInt, string boundName,
                                           float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[SPECIESMLB_BINDING_INDEX]->getSpecies().getName();
          
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void LinkerBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (cc->getCMonomer(bindingSite)->speciesBound(L_BINDING_INDEX)->getN() == 1 ) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getFilament() == c->getFilament()) continue;
            
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites.begin();
                     it != SysParams::Chemistry().bindingSites.end(); it++) {
                
                if (ccn->getCMonomer(*it)->speciesBound(L_BINDING_INDEX)->getN() == 1) {
                    
                    //check distances..
                    auto mp1 = (float)bindingSite / SysParams::Geometry().cylinderIntSize;
                    auto mp2 = (float)*it / SysParams::Geometry().cylinderIntSize;
                    
                    auto x1 = c->getFirstBead()->coordinate;
                    auto x2 = c->getSecondBead()->coordinate;
                    auto x3 = cn->getFirstBead()->coordinate;
                    auto x4 = cn->getSecondBead()->coordinate;
                    
                    auto m1 = midPointCoordinate(x1, x2, mp1);
                    auto m2 = midPointCoordinate(x3, x4, mp2);
                    
                    double dist = twoPointDistance(m1, m2);
                    
                    if(dist > _rMax || dist < _rMin) continue;
                    
                    auto t1 = tuple<CCylinder*, short>(cc, bindingSite);
                    auto t2 = tuple<CCylinder*, short>(ccn, *it);
                    
                    //add in correct order
                    if(c->getID() > cn->getID())
                        _possibleBindings.emplace(t1,t2);
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {
                            
                            _possibleBindings.emplace(t2,t1);
                        }
                        //add in other
                        else {
                            auto m = (LinkerBindingManager*)cn->getCompartment()->
                                      getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);
                            
                            m->_possibleBindings.emplace(t2,t1);
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
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
    
    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
        
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
            if (cc->getCMonomer(*it1)->speciesBound(L_BINDING_INDEX)->getN() == 1 ) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
                    
                    if(cn->getFilament() == c->getFilament()) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites.begin();
                             it2 != SysParams::Chemistry().bindingSites.end(); it2++) {
                        
                        if (ccn->getCMonomer(*it2)->speciesBound(L_BINDING_INDEX)->getN() == 1) {
                            
                            //check distances..
                            auto mp1 = (float)*it1 / SysParams::Geometry().cylinderIntSize;
                            auto mp2 = (float)*it2 / SysParams::Geometry().cylinderIntSize;
                            
                            auto x1 = c->getFirstBead()->coordinate;
                            auto x2 = c->getSecondBead()->coordinate;
                            auto x3 = cn->getFirstBead()->coordinate;
                            auto x4 = cn->getSecondBead()->coordinate;
                            
                            auto m1 = midPointCoordinate(x1, x2, mp1);
                            auto m2 = midPointCoordinate(x3, x4, mp2);
                            
                            double dist = twoPointDistance(m1,m2);
                            
                            if(dist > _rMax || dist < _rMin) continue;
                            
                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);
                            
                            //add in correct order
                            if(c->getID() > cn->getID())
                                _possibleBindings.emplace(t1, t2);
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

MotorBindingManager::MotorBindingManager(ReactionBase* reaction,
                                         Compartment* compartment,
                                         short boundInt, string boundName,
                                         float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
     _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[SPECIESMLB_BINDING_INDEX]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void MotorBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (cc->getCMonomer(bindingSite)->speciesBound(M_BINDING_INDEX)->getN() == 1) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getFilament() == c->getFilament()) continue;
        
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites.begin();
                     it != SysParams::Chemistry().bindingSites.end(); it++) {
                
                if (ccn->getCMonomer(*it)->speciesBound(M_BINDING_INDEX)->getN() == 1) {
                    
                    //check distances..
                    auto mp1 = (float)bindingSite / SysParams::Geometry().cylinderIntSize;
                    auto mp2 = (float)*it / SysParams::Geometry().cylinderIntSize;
                    
                    auto x1 = c->getFirstBead()->coordinate;
                    auto x2 = c->getSecondBead()->coordinate;
                    auto x3 = cn->getFirstBead()->coordinate;
                    auto x4 = cn->getSecondBead()->coordinate;
                    
                    auto m1 = midPointCoordinate(x1, x2, mp1);
                    auto m2 = midPointCoordinate(x3, x4, mp2);
                    
                    double dist = twoPointDistance(m1, m2);
                    
                    if(dist > _rMax || dist < _rMin) continue;
                    
                    auto t1 = tuple<CCylinder*, short>(cc, bindingSite);
                    auto t2 = tuple<CCylinder*, short>(ccn, *it);
                    
                    //add in correct order
                    if(c->getID() > cn->getID())
                        _possibleBindings.emplace(t1,t2);
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {
                            
                            _possibleBindings.emplace(t2,t1);
                        }
                        //add in other
                        else {
                            
                            auto m = (MotorBindingManager*)cn->getCompartment()->
                                      getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);
                            
                            m->_possibleBindings.emplace(t2,t1);
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
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
    
    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
        
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
            if (cc->getCMonomer(*it1)->speciesBound(M_BINDING_INDEX)->getN() == 1) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
                    
                    if(cn->getFilament() == c->getFilament()) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites.begin();
                             it2 != SysParams::Chemistry().bindingSites.end(); it2++) {
                        
                        if (ccn->getCMonomer(*it2)->speciesBound(M_BINDING_INDEX)->getN() == 1) {
                            
                            //check distances..
                            auto mp1 = (float)*it1 / SysParams::Geometry().cylinderIntSize;
                            auto mp2 = (float)*it2 / SysParams::Geometry().cylinderIntSize;
                            
                            auto x1 = c->getFirstBead()->coordinate;
                            auto x2 = c->getSecondBead()->coordinate;
                            auto x3 = cn->getFirstBead()->coordinate;
                            auto x4 = cn->getSecondBead()->coordinate;
                            
                            auto m1 = midPointCoordinate(x1, x2, mp1);
                            auto m2 = midPointCoordinate(x3, x4, mp2);
                            
                            double dist = twoPointDistance(m1, m2);
                            
                            if(dist > _rMax || dist < _rMin) continue;
                            
                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);
                            
                            //add in correct order
                            if(c->getID() > cn->getID())
                                _possibleBindings.emplace(t1,t2);
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

