
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

#include "Cylinder.h"

mt19937* FilamentBindingManager::_eng = 0;

void BranchingManager::updatePossibleBindings(CCylinder* cc, short bindingSite) {
    
    //remove all tuples which have this ccylinder
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));

    //now re add valid binding site
    if (cc->getCMonomer(bindingSite)->activeSpeciesBound() == BOUND_EMPTY)
            _possibleBindings.insert(tuple<CCylinder*, short>(cc, bindingSite));
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
    _possibleBindings.clear();
    
    for(auto &c : _compartment->getCylinders()) {
    
        auto cc = c->getCCylinder();
        
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites.begin();
                 it != SysParams::Chemistry().bindingSites.end(); it++)
            
            if (cc->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY)
                
                _possibleBindings.insert(tuple<CCylinder*, short>(cc, *it));
    }
}

//void LinkerBindingManager::updatePossibleBindings(CCylinder* cc, short bindingSite) {
//    
//    vector<tuple<CCylinder*, short>> toRemove;
//    
//    //remove all tuples which have this ccylinder (only as key)
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//        
//        if (get<0>(it->first) == cc)
//            _possibleBindings.erase(it++);
//        
//        else if(get<0>(it->second) == cc) {
//            
//            toRemove.push_back
//            
//        }
//        
//        else ++it;
//    }
//    //now re add valid binding sites
//    for(auto it = SysParams::Chemistry().bindingSites.begin();
//             it != SysParams::Chemistry().bindingSites.end(); it++) {
//        
//        if (cc->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY) {
//            
//            //loop through neighbors
//            //now re add valid based on CCNL
//            for (auto cn : _neighborList->getNeighbors(cc->getCylinder())) {
//                
//                auto ccn = cn->getCCylinder();
//                
//                for(auto it2 = SysParams::Chemistry().bindingSites.begin();
//                         it2 != SysParams::Chemistry().bindingSites.end(); it2++)
//                    
//                    if (ccn->getCMonomer(*it2)->activeSpeciesBound() == BOUND_EMPTY)
//                        
//                        _possibleBindings.insert(tuple<CCylinder*, short>(cc, *it),
//                                                 tuple<CCylinder*, short>(ccn, *it2));
//            }
//        }
//    }
//}
//
//void LinkerBindingManager::replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) {
//    
//    unordered_map<tuple<CCylinder*, short>, tuple<CCylinder*, short>> toAdd;
//    
//    //remove all tuples which have this ccylinder
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//        
//        if (get<0>(it->first) == oldcc) {
//            
//            //add new tuple
//            toAdd[tuple<CCylinder*, short>(newcc, get<1>(it->first))] = it->second;
//            //erase
//            _possibleBindings.erase(it++);
//        }
//        else if (get<0>(it->second) == oldcc) {
//            
//            //add new tuple
//            toAdd[it->second] = tuple<CCylinder*, short>(newcc, get<1>(it->first));
//            //erase
//            _possibleBindings.erase(it++);
//        }
//        else ++it;
//    }
//    
//    //add all to set
//    for(auto it : toAdd) _possibleBindings[it.first] = it.second;
//}
//
//
//void MotorBindingManager::updatePossibleBindings(CCylinder* cc) {
//    
//    //remove all tuples which have this ccylinder (only as key)
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//        
//        if (get<0>(it->first) == cc)
//            _possibleBindings.erase(it++);
//        
//        else ++it;
//    }
//    //now re add valid binding sites
//    for(auto it = SysParams::Chemistry().bindingSites.begin();
//             it != SysParams::Chemistry().bindingSites.end(); it++) {
//        
//        if (cc->getCMonomer(*it)->activeSpeciesBound() == BOUND_EMPTY) {
//            
//            //loop through neighbors
//            //now re add valid based on CCNL
//            for (auto cn : _neighborList->getNeighbors(cc->getCylinder())) {
//                
//                auto ccn = cn->getCCylinder();
//                
//                for(auto it2 = SysParams::Chemistry().bindingSites.begin();
//                         it2 != SysParams::Chemistry().bindingSites.end(); it2++)
//                    
//                    if (ccn->getCMonomer(*it2)->activeSpeciesBound() == BOUND_EMPTY)
//                        
//                        _possibleBindings.insert(tuple<CCylinder*, short>(cc, *it),
//                                                 tuple<CCylinder*, short>(ccn, *it2));
//            }
//        }
//    }
//}
//
//void MotorBindingManager::replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) {
//    
//    unordered_map<tuple<CCylinder*, short>, tuple<CCylinder*, short>> toAdd;
//    
//    //remove all tuples which have this ccylinder
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//        
//        if (get<0>(it->first) == oldcc) {
//            
//            //add new tuple
//            toAdd[tuple<CCylinder*, short>(newcc, get<1>(it->first))] = it->second;
//            //erase
//            _possibleBindings.erase(it++);
//        }
//        else if (get<0>(it->second) == oldcc) {
//            
//            //add new tuple
//            toAdd[it->second] = tuple<CCylinder*, short>(newcc, get<1>(it->first));
//            //erase
//            _possibleBindings.erase(it++);
//        }
//        else ++it;
//    }
//    //add all to set
//    for(auto it : toAdd) _possibleBindings[it.first] = it.second;
//}
