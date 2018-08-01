
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BindingManager.h"

#include "Compartment.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "MotorGhost.h"

#include "SubSystem.h"
#include "Boundary.h"
#include "CompartmentGrid.h"

#include "ChemCallbacks.h"
#include "MathFunctions.h"
#include "GController.h"
#include "SysParams.h"

using namespace mathfunc;

//BRANCHER

BranchingManager::BranchingManager(ReactionBase* reaction,
                                   Compartment* compartment,
                                   short boundInt, string boundName,
                                   short filamentType,
                                   NucleationZoneType zone, double nucleationDistance)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
      _nucleationZone(zone), _nucleationDistance(nucleationDistance) {
    
    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[B_RXN_INDEX]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
}

void BranchingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
    bool inZone = true;
    //see if in nucleation zone
    if(_nucleationZone != NucleationZoneType::ALL) {
        
        auto mp = (float)bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
        
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
            //Qin, add SIDEBOUNDARY that check the distance to the side of cylinder
            else if(_nucleationZone == NucleationZoneType::SIDEBOUNDARY){
                if(_subSystem->getBoundary()->sidedistance(coord) < _nucleationDistance)
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
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(), 1.0) && inZone) {
        
        auto t = tuple<CCylinder*, short>(cc, bindingSite);
        _possibleBindings.insert(t);
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void BranchingManager::addPossibleBindings(CCylinder* cc) {
    
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        addPossibleBindings(cc, *bit);
}

void BranchingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
    //remove tuple which has this ccylinder
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}


void BranchingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}


void BranchingManager::updateAllPossibleBindings() {
    
    //clear all
    _possibleBindings.clear();
    
    for(auto &c : _compartment->getCylinders()) {
    
        if(c->getType() != _filamentType) continue;
        
        auto cc = c->getCCylinder();
        
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
            
            bool inZone = true;
            //see if in nucleation zone
            if(_nucleationZone != NucleationZoneType::ALL) {
                
                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
                
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
                    //Qin, add SIDEBOUNDARY that check the distance to the side of cylinder
                    else if(_nucleationZone == NucleationZoneType::SIDEBOUNDARY){
                        if(_subSystem->getBoundary()->sidedistance(coord) < _nucleationDistance){
                            inZone = true;
                            //cout << "x= " << coord[1] << "y= " << coord[2] << endl;
                        }
                        
                        
                        else
                            inZone = false;
                    }
                    else inZone = true;
                }
                else
                    inZone = false;
            }
            if (areEqual(cc->getCMonomer(*it)->speciesBound(
                SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(), 1.0) && inZone) {
                
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindings.insert(t);
            }
        }
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

bool BranchingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {
        
        CCylinder* cc = get<0>(*it);
        Cylinder* c   = cc->getCylinder();
        
        short bindingSite = get<1>(*it);
        
        bool flag = true;
    
        //check site empty
        if(!areEqual(cc->getCMonomer(bindingSite)->speciesBound(
           SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(), 1.0))
            flag = false;

        if(!flag) {
            cout << "Binding site in branching manager is inconsistent. " << endl;
            cout << "Binding site = " << bindingSite << endl;
            
            cout << "Cylinder info ..." << endl;
            c->printSelf();
            
            return false;
        }
    }
    return true;
}

//LINKER

LinkerBindingManager::LinkerBindingManager(ReactionBase* reaction,
                                           Compartment* compartment,
                                           short boundInt, string boundName,
                                           short filamentType,
                                           float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
      _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[ML_RXN_INDEX]->getSpecies().getName();
          
    _bindingSpecies = _compartment->findSpeciesByName(name);
}


void LinkerBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
#ifdef DEBUGCONSTANTSEED
    struct Orderset
    {
        bool operator()(Cylinder* lhs, Cylinder* rhs) const  {
            return lhs->getID() < rhs->getID();
        }
    };
#endif
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        vector<Cylinder*> nList = _neighborLists[_nlIndex]->getNeighbors
                (cc->getCylinder());
#ifdef DEBUGCONSTANTSEED
        sort(nList.begin(),nList.end(),Orderset());
#endif
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : nList) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
            
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                
                if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {
                    
                    //check distances..
                    auto mp1 = (float)bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
                    auto mp2 = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
                    
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
                    if(c->getID() > cn->getID()) {
#ifdef DEBUGCONSTANTSEED
                        appendpossibleBindings(t1,t2);
//                        _possibleBindings.emplace(t1, t2);
#else
                        _possibleBindings.emplace(t1, t2);
#endif
                    }
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {

#ifdef DEBUGCONSTANTSEED
                            appendpossibleBindings(t1,t2);
//                            _possibleBindings.emplace(t1, t2);
#else
                            _possibleBindings.emplace(t1, t2);
#endif
                        }
                        //add in other
                        else {
                            auto m = (LinkerBindingManager*)cn->getCompartment()->
                                      getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);
#ifdef DEBUGCONSTANTSEED
                            m->appendpossibleBindings(t2,t1);
//                            m->_possibleBindings.emplace(t2,t1);
#else
                            m->_possibleBindings.emplace(t2,t1);
#endif
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
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
    
        addPossibleBindings(cc, *bit);
}


void LinkerBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

#ifdef DEBUGCONSTANTSEED
    erasepossibleBindings(cc,bindingSite);
    //remove all tuples which have this ccylinder as key
#else
    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
#endif
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
//    std::cout<<"removed rxns "<<abs(oldN-newN)<<"(Old "<<oldN<<" NewN "<<newN<<")"<<endl;
    updateBindingReaction(oldN, newN);
    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
        
        if(cn->getType() != _filamentType) continue;
        
        if(cn->getCompartment() != _compartment) {
            
            auto m = (LinkerBindingManager*)cn->getCompartment()->
                      getFilamentBindingManagers()[_mIndex].get();
            
            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
    //remove, update affected
    for(auto m : affectedManagers) {
#ifdef DEBUGCONSTANTSEED
        m->erasepossibleBindings(cc,bindingSite);
#else
        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {
            
            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindings.erase(it++);
            else ++it;
        }
#endif
        
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        std::cout<<"removed affected mgr rxns "<<abs(oldN-newN)<<"(Old "<<oldN<<" NewN "
                ""<<newN<<endl;
        m->updateBindingReaction(oldNOther, newNOther);
    }
}

void LinkerBindingManager::removePossibleBindings(CCylinder* cc) {
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void LinkerBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
#ifdef DEBUGCONSTANTSEED
    struct Orderset
    {
        bool operator()(Cylinder* lhs, Cylinder* rhs) const  {
            return lhs->getID() < rhs->getID();
        }
    };
    set<Cylinder*, Orderset> _cylinderssorted; ///< Set of cylinders that are in this
    for(auto c : _compartment->getCylinders()) {
        _cylinderssorted.insert(c);
    }
    
    for(auto c : _cylinderssorted)
#else
        for(auto c : _compartment->getCylinders())
#endif
        {
    
        if(c->getType() != _filamentType) continue;
        
        auto cc = c->getCCylinder();
    
        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        
            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0) ) {

                vector<Cylinder*> nList = _neighborLists[_nlIndex]->getNeighbors
                        (cc->getCylinder());
#ifdef DEBUGCONSTANTSEED
                sort(nList.begin(),nList.end(),Orderset());
#endif
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : nList) {
                    
                    if(cn->getParent() == c->getParent()) continue;
                    if(cn->getType() != _filamentType) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                             it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
                        
                        if (areEqual(ccn->getCMonomer(*it2)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {
                            
                            //check distances..
                            auto mp1 = (float)*it1 / SysParams::Geometry().cylinderNumMon[_filamentType];
                            auto mp2 = (float)*it2 / SysParams::Geometry().cylinderNumMon[_filamentType];
                            
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
                            if(c->getID() > cn->getID()) {
#ifdef DEBUGCONSTANTSEED
                                appendpossibleBindings(t1,t2);
#else
                                _possibleBindings.emplace(t1, t2);
#endif
                            }
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

bool LinkerBindingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {
#ifdef DEBUGCONSTANTSEED
        CCylinder* cc1 = get<0>((*it)[0]);
        CCylinder* cc2 = get<0>((*it)[1]);

        short bindingSite1 = get<1>((*it)[0]);
        short bindingSite2 = get<1>((*it)[1]);
#else
        CCylinder* cc1 = get<0>(it->first);

        CCylinder* cc2 = get<0>(it->second);

        short bindingSite1 = get<1>(it->first);
        short bindingSite2 = get<1>(it->second);
#endif
        Cylinder*  c1  = cc1->getCylinder();
        Cylinder*  c2  = cc2->getCylinder();
        bool flag = true;
        
        //check site empty
        if(!areEqual(cc1->getCMonomer(bindingSite1)->speciesBound(
           SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0) ||
           
           !areEqual(cc2->getCMonomer(bindingSite2)->speciesBound(
           SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0))

            flag = false;
        
        if(!flag) {
            cout << "Binding site in linker manager is inconsistent. " << endl;
            cout << "Binding site for cylinder 1 = " << bindingSite1 << endl;
            cout << "Binding site for cylinder 2 = " << bindingSite2 << endl;
            
            cout << "Cylinder info ..." << endl;
            c1->printSelf();
            c2->printSelf();
            
            return false;
        }
    }
    return true;
}

/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> LinkerBindingManager::chooseBindingSites() {

  assert((_possibleBindings.size() != 0)
	            && "Major bug: Linker binding manager should not have zero binding \
           sites when called to choose a binding site.");
    std::cout<<"Choose binding sites"<<endl;
  int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);

    auto it = _possibleBindings.begin();
//
    advance(it, randomIndex);

  auto xxx = _compartment->coordinates();
    std::cout<<"Compartment coords "<<xxx[0]<<" "<<xxx[1]<<" "<<xxx[2]<<endl;
#ifdef DEBUGCONSTANTSEED
    return (*it);
#else
  return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}

#ifdef DEBUGCONSTANTSEED
void LinkerBindingManager::erasepossibleBindings(CCylinder* cc, short bindingSite) {
    tuple<CCylinder*, short> bindingtoremove = make_tuple(cc, bindingSite);
//    std::cout<<"erasing cyl "<<cc->getCylinder()->getID()<<" bs "<<bindingSite<<endl;
    int counter = 0;
    for (auto p = _possibleBindings.begin(); p != _possibleBindings.end(); p++) {
        auto binding1 = (*p)[0];
        auto binding2 = (*p)[1];
        if (bindingtoremove == binding1 || bindingtoremove == binding2)
        {
//            auto cyl1 = get<0>(binding1)->getCylinder()->getID();
//            auto cyl2 = get<0>(binding2)->getCylinder()->getID();
//            auto bs1 = get<1>(binding1);
//            auto bs2 = get<1>(binding2);
//            std::cout<<"Removing pair Cyl "<<cyl1<<" bs "<<bs1<<" Cyl "<<cyl2<<" bs "
//                    ""<<bs2<<"  ID "<<counter<<endl;
            _possibleBindings.erase(p);
        }
        counter++;
    }
}
#endif
void LinkerBindingManager::appendpossibleBindings(tuple<CCylinder*, short> t1,
                                              tuple<CCylinder*,
        short> t2){
    double oldN=numBindingSites();
#ifdef DEBUGCONSTANTSEED
    vector<tuple<CCylinder*, short>> a = {t1,t2};
    bool status = true;
    for(auto p=_possibleBindings.begin();p!=_possibleBindings.end();p++) {
        auto binding1 = (*p)[0];
        auto binding2 = (*p)[1];
        if (t1 == binding1 && t2 == binding2){
            status = false;
            break;
        }
    }
    if(status)
    {  _possibleBindings.push_back(a);
    }
#else
    _possibleBindings.emplace(t1,t2);
#endif

    double newN=numBindingSites();
    updateBindingReaction(oldN,newN);
}

//MOTOR
MotorBindingManager::MotorBindingManager(ReactionBase* reaction,
                                         Compartment* compartment,
                                         short boundInt, string boundName,
                                         short filamentType,
                                         float rMax, float rMin)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
     _rMin(rMin), _rMax(rMax) {
    
    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[ML_RXN_INDEX]->getSpecies().getName();
    
    _bindingSpecies = _compartment->findSpeciesByName(name);
         
    //initialize ID's based on number of species in compartment
//    int numSpecies = rs[ML_RXN_INDEX + 1]->getSpecies().getN();

    //DEPRECATED AS OF 9/22/16
//    for(int i = 0; i < numSpecies; i++)
//        _unboundIDs.push_back(MotorGhost::_motorGhosts.getID());
         
         
    //attach an rspecies callback to this species
    Species* sd = &(rs[ML_RXN_INDEX + 1]->getSpecies());
         
    UpdateMotorIDCallback mcallback(boundInt);
    ConnectionBlock rcb(sd->connect(mcallback,false));
         
}


void MotorBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
#ifdef DEBUGCONSTANTSEED
    struct Orderset
    {
        bool operator()(Cylinder* lhs, Cylinder* rhs) const  {
            return lhs->getID() < rhs->getID();
        }
    };
#endif
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        vector<Cylinder*> nList = _neighborLists[_nlIndex]->getNeighbors
                (cc->getCylinder());
#ifdef DEBUGCONSTANTSEED
        sort(nList.begin(),nList.end(),Orderset());
#endif
        for (auto cn : nList) {
            
            Cylinder* c = cc->getCylinder();
            
            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
        
            auto ccn = cn->getCCylinder();
            
            for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                
                if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
                    
                    //check distances..
                    auto mp1 = (float)bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
                    auto mp2 = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
                    
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
                    {
#ifdef DEBUGCONSTANTSEED
                        appendpossibleBindings(t1,t2);
//                        _possibleBindings.emplace(t1, t2);
#else
                        _possibleBindings.emplace(t1, t2);
#endif
                    }
                    
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {

#ifdef DEBUGCONSTANTSEED
                            appendpossibleBindings(t1,t2);
//                            _possibleBindings.emplace(t1, t2);
#else
                            _possibleBindings.emplace(t1, t2);
#endif
                        }
                        //add in other
                        else {
                            
                            auto m = (MotorBindingManager*)cn->getCompartment()->
                                      getFilamentBindingManagers()[_mIndex].get();
                            
                            affectedManagers.push_back(m);

#ifdef DEBUGCONSTANTSEED
                            m->appendpossibleBindings(t2,t1);
//                            m->_possibleBindings.emplace(t2,t1);
#else
                            m->_possibleBindings.emplace(t2,t1);
#endif
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

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        addPossibleBindings(cc, *bit);
}


void MotorBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
#ifdef DEBUGCONSTANTSEED
    erasepossibleBindings(cc,bindingSite);
    //remove all tuples which have this ccylinder as key
#else
    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);
    
    //remove all tuples which have this as value
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
        
        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindings.erase(it++);
        
        else ++it;
    }
#endif
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
    
    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
        
        if(cn->getType() != _filamentType) continue;
        
        if(cn->getCompartment() != _compartment) {
            
            auto m = (MotorBindingManager*)cn->getCompartment()->
                      getFilamentBindingManagers()[_mIndex].get();
            
            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
    
    //remove, update affected
    for(auto m : affectedManagers) {
#ifdef DEBUGCONSTANTSEED
        m->erasepossibleBindings(cc,bindingSite);
#else
        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindings.erase(it++);
            
            else ++it;
        }
#endif
        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();
        
        m->updateBindingReaction(oldNOther, newNOther);
    }
}

void MotorBindingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void MotorBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
#ifdef DEBUGCONSTANTSEED
    struct Orderset
    {
        bool operator()(Cylinder* lhs, Cylinder* rhs) const  {
            return lhs->getID() < rhs->getID();
        }
    };
    set<Cylinder*, Orderset> _cylinderssorted; ///< Set of cylinders that are in this
    for(auto c : _compartment->getCylinders()) {
        _cylinderssorted.insert(c);
    }
    
    for(auto c : _cylinderssorted)
#else
        for(auto c : _compartment->getCylinders())
#endif
        {

        if(c->getType() != _filamentType) continue;
        
        auto cc = c->getCCylinder();
        
        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
            
            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {

#ifdef DEBUGCONSTANTSEED
                vector<Cylinder*> nList = _neighborLists[_nlIndex]->getNeighbors
                        (cc->getCylinder());
                sort(nList.begin(),nList.end(),Orderset());
                for (auto cn : nList)
#endif
for(auto cn:_neighborLists[_nlIndex]->getNeighbors
        (cc->getCylinder()))
                //loop through neighbors
                //now re add valid based on CCNL
                 {
                    
                    if(cn->getParent() == c->getParent()) continue;
                    if(cn->getType() != _filamentType) continue;
                    
                    auto ccn = cn->getCCylinder();
                    
                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                             it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
                        
                        if (areEqual(ccn->getCMonomer(*it2)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
                            
                            //check distances..
                            auto mp1 = (float)*it1 / SysParams::Geometry().cylinderNumMon[_filamentType];
                            auto mp2 = (float)*it2 / SysParams::Geometry().cylinderNumMon[_filamentType];
                            
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
                            if(c->getID() > cn->getID()){
#ifdef DEBUGCONSTANTSEED
                                appendpossibleBindings(t1,t2);
#else
                                _possibleBindings.emplace(t1, t2);
#endif
                            }
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

bool MotorBindingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {

#ifdef DEBUGCONSTANTSEED
        CCylinder* cc1 = get<0>((*it)[0]);
        CCylinder* cc2 = get<0>((*it)[1]);

        short bindingSite1 = get<1>((*it)[0]);
        short bindingSite2 = get<1>((*it)[1]);
//        CCylinder* cc1 = get<0>(it->first);
//
//        CCylinder* cc2 = get<0>(it->second);
//
//        short bindingSite1 = get<1>(it->first);
//        short bindingSite2 = get<1>(it->second);
#else
        CCylinder* cc1 = get<0>(it->first);

        CCylinder* cc2 = get<0>(it->second);

        short bindingSite1 = get<1>(it->first);
        short bindingSite2 = get<1>(it->second);
#endif
        Cylinder*  c1  = cc1->getCylinder();
        Cylinder*  c2  = cc2->getCylinder();
        bool flag = true;
        
        //check site empty
        if(!areEqual(cc1->getCMonomer(bindingSite1)->speciesBound(
           SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0) ||
           
           !areEqual(cc2->getCMonomer(bindingSite2)->speciesBound(
           SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0))
            
            flag = false;
        
        if(!flag) {
            cout << "Binding site in motor manager is inconsistent. " << endl;
            cout << "Binding site for cylinder 1 = " << bindingSite1 << endl;
            cout << "Binding site for cylinder 2 = " << bindingSite2 << endl;
            
            cout << "Cylinder info ..." << endl;
            c1->printSelf();
            c2->printSelf();
            
            //check if in neighbor list
            auto nlist = _neighborLists[_nlIndex]->getNeighbors(c1);
            
            if(find(nlist.begin(), nlist.end(), c2) == nlist.end()) {
                
                
                cout << "Not in neighbor list 1" << endl;
                
            }
            
            nlist = _neighborLists[_nlIndex]->getNeighbors(c2);
            
            if(find(nlist.begin(), nlist.end(), c1) == nlist.end()) {
                
                
                cout << "Not in neighbor list 2" << endl;
                
            }
            
            
            return false;
        }
    }
    return true;
}

/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> MotorBindingManager::chooseBindingSites() {

    assert((_possibleBindings.size() != 0)
           && "Major bug: Motor binding manager should not have zero binding \
           sites when called to choose a binding site.");
    std::cout<<"Choose binding sites"<<endl;
    int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
    auto it = _possibleBindings.begin();

    advance(it, randomIndex);
    auto xxx = _compartment->coordinates();
    std::cout<<"Compartment coords "<<xxx[0]<<" "<<xxx[1]<<" "<<xxx[2]<<endl;
#ifdef DEBUGCONSTANTSEED
    return (*it);
#else
    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}
#ifdef DEBUGCONSTANTSEED
void MotorBindingManager::erasepossibleBindings(CCylinder* cc, short bindingSite) {
    tuple<CCylinder*, short> bindingtoremove = make_tuple(cc, bindingSite);
//    std::cout<<"erasing cyl "<<cc->getCylinder()->getID()<<" bs "<<bindingSite<<endl;
    int counter = 0;
    for (auto p = _possibleBindings.begin(); p != _possibleBindings.end(); p++) {
        auto binding1 = (*p)[0];
        auto binding2 = (*p)[1];
        if (bindingtoremove == binding1 || bindingtoremove == binding2)
        {
//            auto cyl1 = get<0>(binding1)->getCylinder()->getID();
//            auto cyl2 = get<0>(binding2)->getCylinder()->getID();
//            auto bs1 = get<1>(binding1);
//            auto bs2 = get<1>(binding2);
//            auto x = _compartment->coordinates();
//            std::cout<<"Removing pair Cyl "<<cyl1<<" bs "<<bs1<<" Cyl "<<cyl2<<" bs "
//                    ""<<bs2<<"  ID "<<counter<<" in compartment "<<x[0]<<" "<<x[1]<<" "
//                    ""<<x[2]<<endl;
            _possibleBindings.erase(p);
        }
        counter++;
    }
}
#endif

void MotorBindingManager::appendpossibleBindings(tuple<CCylinder*, short> t1,
                                              tuple<CCylinder*, short> t2){
    double oldN=numBindingSites();
#ifdef DEBUGCONSTANTSEED
    vector<tuple<CCylinder*, short>> a = {t1,t2};
    bool status = true;
    for(auto p=_possibleBindings.begin();p!=_possibleBindings.end();p++) {
        auto binding1 = (*p)[0];
        auto binding2 = (*p)[1];
        if (t1 == binding1 && t2 == binding2){
            status = false;
            break;
        }
    }
    if(status)
    {  _possibleBindings.push_back(a);
//        auto c1 = (get<0>(t1))->getCylinder()->coordinate;
//        auto c2 = (get<0>(t2))->getCylinder()->coordinate;
//        auto c1ID = (get<0>(t1))->getCylinder()->getID();
//        auto c2ID = (get<0>(t2))->getCylinder()->getID();
//        auto bs1 = get<1>(t1);
//        auto bs2 = get<1>(t2);
//        auto x = _compartment->coordinates();
//        std::cout<<"Added Cyl "<<c1ID<<" bs "<<bs1<<" Cyl "<<c2ID<<" bs "<<bs2<<" "
//                "to position "<<_possibleBindings.size()<<" in compartment "<<x[0]<<" "<<x[1]<<" "
//                         ""<<x[2]<<"with coords "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "
//                         ""<<c2[1]<<" "<<c2[1]<<" "<<c2[2]<<endl;
    }
//    _possibleBindings.emplace(t1,t2);
//            auto c1 = (get<0>(t1))->getCylinder()->coordinate;
//        auto c2 = (get<0>(t2))->getCylinder()->coordinate;
//        auto c1ID = (get<0>(t1))->getCylinder()->getID();
//        auto c2ID = (get<0>(t2))->getCylinder()->getID();
//        auto bs1 = get<1>(t1);
//        auto bs2 = get<1>(t2);
//        auto x = _compartment->coordinates();
//        std::cout<<"Added Cyl "<<c1ID<<" bs "<<bs1<<" Cyl "<<c2ID<<" bs "<<bs2<<" "
//                "to position "<<_possibleBindings.size()<<" in compartment "<<x[0]<<" "<<x[1]<<" "
//                         ""<<x[2]<<"with coords "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "
//                         ""<<c2[1]<<" "<<c2[1]<<" "<<c2[2]<<endl;
#else
    _possibleBindings.emplace(t1,t2);
#endif

    double newN=numBindingSites();
    updateBindingReaction(oldN,newN);
}

SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CylinderCylinderNL*> LinkerBindingManager::_neighborLists;
vector<CylinderCylinderNL*> MotorBindingManager::_neighborLists;

