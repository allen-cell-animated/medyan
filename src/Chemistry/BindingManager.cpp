
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
#include "../Structure/CaMKIICylinder.h"

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

//CAMKII Binding

CaMKIIBindingManager::CaMKIIBindingManager(ReactionBase* reaction,
                                   Compartment* compartment,
                                   short boundInt, string boundName,
                                   short filamentType,
                                   NucleationZoneType zone, double nucleationDistance)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
      _nucleationZone(zone), _nucleationDistance(nucleationDistance) {
    
    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[C1_RXN_INDEX]->getSpecies().getName();

	_bindingSpecies = _compartment->findSpeciesByName(name);
}

void CaMKIIBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

	// Only non-CaMKII cylinders are allowed to run this method.
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
            else inZone = true;
        }
        else
            inZone = false;
    }
    
    //add valid site
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().camkiierBoundIndex[_filamentType])->getN(), 1.0) && inZone) {
        
        auto t = tuple<CCylinder*, short>(cc, bindingSite);
        _possibleBindings.insert(t);
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

void CaMKIIBindingManager::addPossibleBindings(CCylinder* cc) {
    
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        addPossibleBindings(cc, *bit);
}

void CaMKIIBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
    //remove tuple which has this ccylinder
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}


void CaMKIIBindingManager::removePossibleBindings(CCylinder* cc) {
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}


void CaMKIIBindingManager::updateAllPossibleBindings() {
    
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
                    else inZone = true;
                }
                else
                    inZone = false;
            }
            if (areEqual(cc->getCMonomer(*it)->speciesBound(
                SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType])->getN(), 1.0) && inZone) {
                
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindings.insert(t);
            }
        }
    }
    
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    
    updateBindingReaction(oldN, newN);
}

bool CaMKIIBindingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {
        
        CCylinder* cc = get<0>(*it);
        Cylinder* c   = cc->getCylinder();
        
        short bindingSite = get<1>(*it);
        
        bool flag = true;
    
        //check site empty
        if(!areEqual(cc->getCMonomer(bindingSite)->speciesBound(
           SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType])->getN(), 1.0))
            flag = false;

        if(!flag) {
            cout << "Binding site in camkiiing manager is inconsistent. " << endl;
            cout << "Binding site = " << bindingSite << endl;
            
            cout << "Cylinder info ..." << endl;
            c->printSelf();
            
            return false;
        }
    }
    return true;
}

//CAMKII Bundling

CaMKIIBundlingManager::CaMKIIBundlingManager(ReactionBase* reaction,
                                   Compartment* compartment,
                                   short boundInt, string boundName,
                                   short filamentType,
				                   float rMax, float rMin, int maxCoordination)

    : FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
    _rMax(rMax), _rMin(rMin), _maxCoordination(maxCoordination) {

	//find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[C2_RXN_INDEX]->getSpecies().getName();

    _bindingSpecies = _compartment->findSpeciesByName(name);
}

bool CaMKIIBundlingManager::checkMaxCoordinationNum(CCylinder* cc) {
	auto dc = dynamic_cast<CaMKIICylinder*>(cc->getCylinder());
	CaMKIIingPoint *cp;
	if (dc) {
		cp = dc->getCaMKIIPointParent();
		return (cp->getCoordinationNumber() >= _maxCoordination);
	}
	return false;
}

void CaMKIIBundlingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

	// Run the function is input cylinder (cc) is normal or CaMKII cylinder.
	if(cc->getType() != _filamentType && cc->getType() != CAMKII_CYLINDER_FILAMENT_TYPE) return;

	// checkType is the type of the binding cylinder.
	// If the cylinder is a CaMKII cylinder, the binding cylinder is a normal cylinder and vice versa
	short checkType = _filamentType;
	if(cc->getType() == _filamentType) {
		checkType = CAMKII_CYLINDER_FILAMENT_TYPE;
	} else {
		//If it is called with CaMKII only execute this script once (when the bindingSite is 0)
		if (bindingSite!=0) return;
	}

	//if we change other managers copy number
	vector<CaMKIIBundlingManager*> affectedManagers;

	//add valid binding sites
	for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
		Cylinder *c = cc->getCylinder();

		if (cn->getType() != checkType) continue;

		// Assign correct description to cylinders
		Cylinder *camkii_c;
		Cylinder *filament_c;
		if (checkType == CAMKII_CYLINDER_FILAMENT_TYPE) {
			camkii_c = cn;
			filament_c = c;
		} else {
			camkii_c = c;
			filament_c = cn;
		}
		CCylinder *camkii_cc = camkii_c->getCCylinder();
		CCylinder *filament_cc = filament_c->getCCylinder();

		//Check that the CaMKII coordination number is less than the maximum
		if (checkMaxCoordinationNum(camkii_cc)) continue;

		if (checkType == CAMKII_CYLINDER_FILAMENT_TYPE) {
			// It was called with a regular cylinder
			// If the binding site of this cylinder is not free, then continue
			if (!areEqual(filament_cc->getCMonomer(bindingSite)->speciesBound(
					SysParams::Chemistry().camkiierBundlingBoundIndex[_filamentType])->getN(), 1.0))
				continue;
			auto mp1 = (float) bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
			auto x1 = camkii_c->getFirstBead()->coordinate;
			auto x3 = filament_c->getFirstBead()->coordinate;
			auto x4 = filament_c->getSecondBead()->coordinate;
			auto m2 = midPointCoordinate(x3, x4, mp1);
			double dist = twoPointDistance(x1, m2);
			if (dist > _rMax || dist < _rMin) continue;
			auto t1 = tuple<CCylinder *, short>(camkii_cc, 0);
			auto t2 = tuple<CCylinder *, short>(filament_cc, bindingSite);

			// We should update the possible binding of the compartment of the CaMKII.
			auto m = (CaMKIIBundlingManager*)camkii_c->getCompartment()->
					getFilamentBindingManagers()[_mIndex].get();
			m->_possibleBindings.emplace(t1, t2);
		} else {
			// It was called with a CaMKII cylinder, iterate over all possible binding sites of the new cylinder
			for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
				 it != SysParams::Chemistry().bindingSites[_filamentType].end();
				 it++) {
				//If the binding site of the neighbor is not free, then continue
				if (!areEqual(filament_cc->getCMonomer(*it)->speciesBound(
						SysParams::Chemistry().camkiierBundlingBoundIndex[_filamentType])->getN(), 1.0))
					continue;
				auto mp2 = (float) *it / SysParams::Geometry().cylinderNumMon[_filamentType];
				auto x1 = camkii_c->getFirstBead()->coordinate;
				auto x3 = filament_c->getFirstBead()->coordinate;
				auto x4 = filament_c->getSecondBead()->coordinate;
				auto m2 = midPointCoordinate(x3, x4, mp2);
				double dist = twoPointDistance(x1, m2);
				if (dist > _rMax || dist < _rMin) continue;
				auto t1 = tuple<CCylinder *, short>(camkii_cc, bindingSite);
				auto t2 = tuple<CCylinder *, short>(filament_cc, *it);
				_possibleBindings.emplace(t1, t2);
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

	size_t s1 = _possibleBindings.size();
	updateAllPossibleBindings();
	size_t s2 = _possibleBindings.size();
	if((s2 - s1 > 1) || (s2-s1 < 0) ) {
		cerr << "======MILLAD: ERROR in size of possible size. --- " << __FILE__ << "   " << __LINE__ << "  " << __FUNCTION__ << endl;
	}
}

void CaMKIIBundlingManager::addPossibleBindings(CCylinder* cc) {
	if(cc->getType() == _filamentType) {
		for (auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
				bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
			addPossibleBindings(cc, *bit);
	} else {
		addPossibleBindings(cc, 0);
	}

	checkLengthOfPossibleBinding();
}

void CaMKIIBundlingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

	// Run the function is input cylinder (cc) is normal or CaMKII cylinder.
	if(cc->getType() != _filamentType && cc->getType() != CAMKII_CYLINDER_FILAMENT_TYPE) return;

	// checkType is the type of the binding cylinder.
	// If the cylinder is a CaMKII cylinder, the binding cylinder is a normal cylinder and viceversa
	short checkType = _filamentType;
	if(cc->getType() == _filamentType) {
		checkType = CAMKII_CYLINDER_FILAMENT_TYPE;
	} else {
		//If it is called with CaMKII only execute this script once (when the bindingSite is 0)
		if (bindingSite!=0) return;
	}

	if(cc->getType() == CAMKII_CYLINDER_FILAMENT_TYPE) {
		//Remove all the tuples that have this value as a key
		auto t = tuple<CCylinder*, short>(cc, bindingSite);
		_possibleBindings.erase(t);
	} else {
		//Remove all the tuples that have this as a value
		for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
			if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
				_possibleBindings.erase(it++);
			else ++it;
		}
	}

	//if we change other managers copy number
	vector<CaMKIIBundlingManager*> affectedManagers;

	int oldN = _bindingSpecies->getN();
	int newN = numBindingSites();

	updateBindingReaction(oldN, newN);

	//remove all neighbors which have this binding site pair
	for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

		if (cn->getType() != _filamentType && cn->getType() != CAMKII_CYLINDER_FILAMENT_TYPE) continue;

		if(cn->getCompartment() != _compartment) {

			auto m = (CaMKIIBundlingManager*)cn->getCompartment()->
					getFilamentBindingManagers()[_mIndex].get();

			if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
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


void CaMKIIBundlingManager::removePossibleBindings(CCylinder* cc) {
	if(cc->getType() == _filamentType) {
		for (auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
			 bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
			removePossibleBindings(cc, *bit);
	} else {
		removePossibleBindings(cc, 0);
	}

	checkLengthOfPossibleBinding();
}

void CaMKIIBundlingManager::updateCaMKIIPossibleBindings(CCylinder* cc) {
  if (checkMaxCoordinationNum(cc)) removePossibleBindings(cc);
}


int CaMKIIBundlingManager::checkLengthOfPossibleBindingInternal(unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> &mm1, unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> &mm2) {
	multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> sorted_mm1(mm1.begin(), mm1.end()), sorted_mm2(mm2.begin(), mm2.end());
	multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> intersection;
	std::set_intersection(sorted_mm1.begin(), sorted_mm1.end(),
						  sorted_mm2.begin(), sorted_mm2.end(),
						  std::inserter(intersection, intersection.end()));
	int counter = 0;
	if(intersection.size() < mm1.size()) {
		for(auto x : intersection) {
			CCylinder *ccc = get<0>(x.first);
			CaMKIICylinder *camccc = dynamic_cast<CaMKIICylinder*>(ccc->getCylinder());
			if(camccc != nullptr)
				counter++;
		}
	}
	cerr << "=========== MILLAD: checkLengthOfPossibleBindingInternal " << counter << "    " << intersection.size() << endl;
	return intersection.size();
}

void CaMKIIBundlingManager::updateAllPossibleBindings() {

	_possibleBindings.clear();
	int camkiiFilamentType = CAMKII_CYLINDER_FILAMENT_TYPE;

    //The first cylinder should be a CAMKII cylinder
    for(auto c : _compartment->getCylinders()) {

        if(c->getType() != camkiiFilamentType) continue;

        auto cc = c->getCCylinder();
        CaMKIIingPoint *cp = nullptr;
        if (!dynamic_cast<CaMKIICylinder*>(c)) {
        	cerr << "CaMKII - There is a CaMKII cylinder without CaMKIIPoint." << endl;
        	exit(1);
        }

        // Binding-site index of CaMKII cylinder
        short it1 = 0;

        // skip if parent coordination number isn't between >=1 and <6 (MAX coordination number SysParams::Chemistry().maxcamkii_coord_number)
        // now re add valid binding sites
        if (checkMaxCoordinationNum(cc)) continue;


        // loop through neighbors
        // The neighbors should be the cylinders from the other filaments (obtained from the neighbor list)
        // now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
            auto ccn = cn->getCCylinder();

            for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {

                if (areEqual(ccn->getCMonomer(*it2)->speciesBound(SysParams::Chemistry().camkiierBundlingBoundIndex[_filamentType])->getN(), 1.0)) {
                    //check distances..
                    auto mp1 = (float) it1 / SysParams::Geometry().cylinderNumMon[camkiiFilamentType];
                    auto mp2 = (float)*it2 / SysParams::Geometry().cylinderNumMon[_filamentType];

                    auto x1 = c->getFirstBead()->coordinate;
                    auto x2 = c->getSecondBead()->coordinate;
                    auto x3 = cn->getFirstBead()->coordinate;
                    auto x4 = cn->getSecondBead()->coordinate;
                    auto m1 = midPointCoordinate(x1, x2, mp1);
                    auto m2 = midPointCoordinate(x3, x4, mp2);

                    double dist = twoPointDistance(m1,m2);

                    if(dist > _rMax || dist < _rMin) continue;

                    auto t1 = tuple<CCylinder*, short>(cc, it1);
                    auto t2 = tuple<CCylinder*, short>(ccn, *it2);

                    //add in correct order
                    _possibleBindings.emplace(t1, t2);
                }
            }
        }
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}


bool CaMKIIBundlingManager::isConsistent() {

    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {

        CCylinder* cc1 = get<0>(it->first);
        Cylinder*  c1  = cc1->getCylinder();

        CCylinder* cc2 = get<0>(it->second);
        Cylinder*  c2  = cc2->getCylinder();

        short bindingSite1 = get<1>(it->first);
        short bindingSite2 = get<1>(it->second);

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
    
    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
            
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

		    //auto bc1 = cc->getCylinder()->getBranchingCylinder(); //Yossi
                    //auto bc2 = ccn->getCylinder()->getBranchingCylinder(); //Yossi
		    auto bc3 = cc->getCylinder()->getMotherCylinder(); //Carlos
		    auto bc4 = ccn->getCylinder()->getMotherCylinder(); //Carlos

		    if(bc3 || bc4) {
		 
 		    continue;} //Yossi
		                   

                    
 
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
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
    
        addPossibleBindings(cc, *bit);
}


void LinkerBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
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
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void LinkerBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
    
        if(c->getType() != _filamentType) continue;
        
        auto cc = c->getCCylinder();
    
        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        
            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0) ) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
                    
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
			    
			    //auto bc1 = cc->getCylinder()->getBranchingCylinder(); //Carlos
		            //auto bc2 = ccn->getCylinder()->getBranchingCylinder(); //Carlos
			    auto bc3 = cc->getCylinder()->getMotherCylinder(); //Carlos
		            auto bc4 = ccn->getCylinder()->getMotherCylinder(); //Carlos

		    	    if(bc3 || bc4) {
 		            continue;} //Carlos
		                   
                            
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

bool LinkerBindingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {
        
        CCylinder* cc1 = get<0>(it->first);
        Cylinder*  c1  = cc1->getCylinder();
        
        CCylinder* cc2 = get<0>(it->second);
        Cylinder*  c2  = cc2->getCylinder();
        
        short bindingSite1 = get<1>(it->first);
        short bindingSite2 = get<1>(it->second);
        
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
    int numSpecies = rs[ML_RXN_INDEX + 1]->getSpecies().getN();

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
    
    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;
    
    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
        SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
        
        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
            
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
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        addPossibleBindings(cc, *bit);
}


void MotorBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {
    
    if(cc->getType() != _filamentType) return;
    
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
    
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
             bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        
        removePossibleBindings(cc, *bit);
}

void MotorBindingManager::updateAllPossibleBindings() {
    
    _possibleBindings.clear();
    
    for(auto c : _compartment->getCylinders()) {
        
        if(c->getType() != _filamentType) continue;
        
        auto cc = c->getCCylinder();
        
        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
            
            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
                
                //loop through neighbors
                //now re add valid based on CCNL
                for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {
                    
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

bool MotorBindingManager::isConsistent() {
    
    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++) {
        
        CCylinder* cc1 = get<0>(it->first);
        Cylinder*  c1  = cc1->getCylinder();
        
        CCylinder* cc2 = get<0>(it->second);
        Cylinder*  c2  = cc2->getCylinder();
        
        short bindingSite1 = get<1>(it->first);
        short bindingSite2 = get<1>(it->second);
        
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

SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CylinderCylinderNL*> LinkerBindingManager::_neighborLists;
vector<CylinderCylinderNL*> MotorBindingManager::_neighborLists;
//vector<CaMKIIingPointCylinderNL*> CaMKIIBundlingManager::_neighborLists; //TODO fix the type on NL
vector<CylinderCylinderNL*> CaMKIIBundlingManager::_neighborLists; //TODO fix the type on NL

