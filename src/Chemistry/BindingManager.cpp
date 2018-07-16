
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
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif


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
#ifdef NLORIGINAL
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
#endif

void BranchingManager::addPossibleBindings(CCylinder* cc) {


    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
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
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        removePossibleBindingsstencil(cc, *bit);
#endif
    }
}

#ifdef NLORIGINAL
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
//                output test
//                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
//                auto x1 = cc->getCylinder()->getFirstBead()->coordinate;
//                auto x2 = cc->getCylinder()->getSecondBead()->coordinate;
//
//                auto coord = midPointCoordinate(x1, x2, mp);
//                std::cout<<c->_dcIndex<<" "<<*it<<" "<<_subSystem->getBoundary()->distance(coord)<<endl;
//                end
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindings.insert(t);
            }
        }
    }
//        std::cout<<_possibleBindings.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}
#endif

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

#ifdef NLSTENCILLIST
void BranchingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {

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
        _possibleBindingsstencil.insert(t);
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}
void BranchingManager::updateAllPossibleBindingsstencil() {

    //clear all
    _possibleBindingsstencil.clear();

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
//                output test
//                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
//                auto x1 = cc->getCylinder()->getFirstBead()->coordinate;
//                auto x2 = cc->getCylinder()->getSecondBead()->coordinate;
//
//                auto coord = midPointCoordinate(x1, x2, mp);
//                std::cout<<c->_dcIndex<<" "<<*it<<" "<<_subSystem->getBoundary()->distance(coord)<<endl;
//                end
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindingsstencil.insert(t);
            }
        }
    }
//        std::cout<<_possibleBindings.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}
void BranchingManager::removePossibleBindingsstencil(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)

        removePossibleBindingsstencil(cc, *bit);
}
void BranchingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //remove tuple which has this ccylinder
    _possibleBindingsstencil.erase(tuple<CCylinder*, short>(cc, bindingSite));

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}

void BranchingManager::crosscheck(){
    //cout<<"Branching NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
    //       <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size())
        cout<<"Branching.. The two methods compared do not yield the same number of "
                "binding sites"<<endl;
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;
        bool state = false;
        auto cyl1_o = get<0>(*it1)->getCylinder();
        auto bs1_o =  get<1>(*it1);
//        auto cyl2_o = get<0>(*it1[1])->getCylinder();
//        auto bs2_o =  get<1>(*it1[1]);
        //
        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(*it2)->getCylinder();
            auto bs1_s =  get<1>(*it2);
//            auto cyl2_s = get<0>(*it2[1])->getCylinder();
//            auto bs2_s =  get<1>(*it2[1]);
            sum = 0;
            if(cyl1_o->getID() == cyl1_s->getID() )
                sum++;
//            if(cyl2_o->getID() == cyl2_s->getID() )
//                sum++;
            if(bs1_o == bs1_s )
                sum++;
//            if(bs2_o == bs2_s )
//                sum++;
            if(sum == 2) {
//                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else {
                    cout << "ERROR. Multiple matches for chosen binding site pair in "
                            "STENCILLIST. Check stencillist.." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        }
    std::cout<<"possible bindings size "<<_possibleBindings.size()<<" Total matches "<<
                                                                                     matches<<endl;
    if(_possibleBindings.size() != matches || _possibleBindings.size() !=
                                                      _possibleBindingsstencil.size()){
        cout<<"Branching. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
#endif
#ifdef CUDAACCL_NL
void BranchingManager::assigncudavars() {
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_distance, 2 * sizeof(double)),"cuda data transfer", " "
            "BindingManager.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)),"cuda data transfer", " "
            "BindingManager.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_zone, sizeof(int)),"cuda data transfer", " "
            "BindingManager.cu");
    double dist[2];
    dist[0] = _nucleationDistance;
    dist[1] = GController::getCenter()[2];
    int n[1];
    n[0] = 0;
    int zone[1];
    zone[0] = -1;
    if(_nucleationZone == NucleationZoneType::ALL)
        zone[0] =0;
    else if(_nucleationZone == NucleationZoneType::BOUNDARY)
        zone[0] =1;
    else if(_nucleationZone == NucleationZoneType::TOPBOUNDARY)
        zone[0] =2;
    CUDAcommon::handleerror(cudaMemcpy(gpu_distance, dist, 2 * sizeof(double), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_zone, zone, sizeof(int), cudaMemcpyHostToDevice));
//    delete dist;
}

void BranchingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_distance),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_zone),"cudaFree", "BindingManager");
}


int* BranchingManager::getzoneCUDA(){
    return gpu_zone;
}
int* BranchingManager::getnumpairsCUDA(){
    return gpu_numpairs;
}
#endif

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

#ifdef NLORIGINAL
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
#endif

void LinkerBindingManager::addPossibleBindings(CCylinder* cc) {

//    auto cylcoord = cc->getCylinder()->coordinate;
//    std::cout<<"L Adding possible to Cyl "<<cc->getCylinder()->getID()<<" w coords "
//            ""<<cylcoord[0]<<" "<<cylcoord[1]<<" "<<cylcoord[2]<<endl;

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
}

void LinkerBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

#ifdef DEBUGCONSTANTSEED
    erasepossibleBindings(cc,bindingSite);
    //remove all tuples which have this ccylinder as key
//    auto t = tuple<CCylinder*, short>(cc, bindingSite);
//    _possibleBindings.erase(t);
//
//    //remove all tuples which have this as value
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//
//        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
//            _possibleBindings.erase(it++);
//
//        else ++it;
//    }
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
//        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {
//
//            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
//                m->_possibleBindings.erase(it++);
//            else ++it;
//        }
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

void LinkerBindingManager::removePossibleBindings(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        removePossibleBindingsstencil(cc, *bit);
#endif
    }
}
#ifdef NLORIGINAL
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
    for(auto c :_cylinderssorted)
#else
        for(auto c :_compartment->getCylinders())
#endif
    {
//        auto x = _compartment->coordinates();
//        std::cout<<"L updateall Cyl "<<c->getID()<<" "<<c->coordinate[0]<<" "<<
//                 c->coordinate[1]<<" "<<c->coordinate[2]<<" in comp "<<x[0]<<" "
//                         ""<<x[1]<<" "<<x[2]<<endl;

        if(c->getType() != _filamentType) continue;

        auto cc = c->getCCylinder();

        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {

            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0) ) {

                //loop through neighbors
                //now re add valid based on CCNL
                vector<Cylinder*> nList = _neighborLists[_nlIndex]->getNeighbors
                        (cc->getCylinder());
#ifdef DEBUGCONSTANTSEED
                sort(nList.begin(),nList.end(),Orderset());
#endif
                for (auto cn :nList) {

                    if(cn->getParent() == c->getParent()) continue;
                    if(cn->getType() != _filamentType) continue;

                    auto ccn = cn->getCCylinder();
//                                                std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<endl;

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
//                            std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<" "<<*it1<<" "<<*it2<<" "<<m1[0]<<" "
//                                    ""<<m1[1]<<" "<<m1[2]<<" "<<m2[0]<<" "<<m2[1]<<" "<<m2[2]<<" "<<mp1<<" "<<mp2<<endl;

                            if(dist > _rMax || dist < _rMin) continue;
//                            std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<" "<<*it1<<" "<<*it2<<" "<<m1[0]<<" "
//                                    ""<<m1[1]<<" "<<m1[2]<<" "<<m2[0]<<" "<<m2[1]<<" "<<m2[2]<<" "<<mp1<<" "<<mp2<<endl;
                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);

                            //add in correct order
                            if(c->getID() > cn->getID()) {
#ifdef DEBUGCONSTANTSEED
                                appendpossibleBindings(t1,t2);
//                                std::cout<<"placing L Cyl "<<c->getID()<<" bs "<<*it1<<" "
//                                        "Cyl "<<cn->getID()<<" bs "<<*it2<<endl;
//                                _possibleBindings.emplace(t1, t2);
#else
                                _possibleBindings.emplace(t1, t2);
#endif
                            }
                        }
                    }
//                    std::cout<<_possibleBindings.size()<<endl;
                }
            }
        }
    }
//    auto x = _compartment->coordinates();
//    std::cout<<"L possibile "<<_possibleBindings.size()<<" in comp "<<x[0]<<" "
//            ""<<x[1]<<" "<<x[2]<<endl;
//    std::cout<<"*****"<<endl;
//    for(auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++){
//        auto bspair = *it;
//        auto cyl1 = get<0>(bspair[0])->getCylinder()->getID();
//        auto cyl2 = get<0>(bspair[1])->getCylinder()->getID();
//        std::cout<<"Cyl "<<cyl1<<"bs "<<get<1>(bspair[0])<<" Cyl "<<cyl2<<" bs "
//                ""<<get<1>(bspair[1])<<endl;
//    }
//    for(auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++){
//        auto cyl1 = get<0>(it->first)->getCylinder();
//        auto cyl2 = get<0>(it->second)->getCylinder();
//        std::cout<<"Cyl "<<cyl1->getID()<<"bs "<<get<1>(it->first)<<" Cyl "<<cyl2->getID
//                ()<<" bs "<<get<1>(it->second)<<endl;
//    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
#ifdef CUDAACCL_NL

#endif
    updateBindingReaction(oldN, newN);
}
#endif

bool LinkerBindingManager::isConsistent() {

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

#ifdef NLSTENCILLIST
void LinkerBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {

        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

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
                    if(c->getID() > cn->getID())
                        _possibleBindingsstencil.emplace(t1,t2);
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {

                            _possibleBindingsstencil.emplace(t2,t1);
                        }
                            //add in other
                        else {
                            auto m = (LinkerBindingManager*)cn->getCompartment()->
                                    getFilamentBindingManagers()[_mIndex].get();

                            affectedManagers.push_back(m);

                            m->_possibleBindingsstencil.emplace(t2,t1);
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
void LinkerBindingManager::updateAllPossibleBindingsstencil() {

    _possibleBindingsstencil.clear();

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
//                                                std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<endl;

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
//                            std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<" "<<*it1<<" "<<*it2<<" "<<m1[0]<<" "
//                                    ""<<m1[1]<<" "<<m1[2]<<" "<<m2[0]<<" "<<m2[1]<<" "<<m2[2]<<" "<<mp1<<" "<<mp2<<endl;

                            if(dist > _rMax || dist < _rMin) continue;
//                            std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<" "<<*it1<<" "<<*it2<<" "<<m1[0]<<" "
//                                    ""<<m1[1]<<" "<<m1[2]<<" "<<m2[0]<<" "<<m2[1]<<" "<<m2[2]<<" "<<mp1<<" "<<mp2<<endl;
                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);

                            //add in correct order
                            if(c->getID() > cn->getID())
                                _possibleBindingsstencil.emplace(t1, t2);
                        }
                    }
//                    std::cout<<_possibleBindings.size()<<endl;
                }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
#ifdef CUDAACCL_NL

#endif
    updateBindingReaction(oldN, newN);
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)

        removePossibleBindingsstencil(cc, *bit);
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencil.erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencil.begin(); it != _possibleBindingsstencil.end();
            ) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindingsstencil.erase(it++);

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

        for (auto it = m->_possibleBindingsstencil.begin(); it !=
                m->_possibleBindingsstencil.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindingsstencil.erase(it++);

            else ++it;
        }

        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();

        m->updateBindingReaction(oldNOther, newNOther);
    }
}

void LinkerBindingManager::crosscheck(){
    //cout<<"Branching NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
    //    <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size())
        cout<<"Linker.. The two methods compared do not yield the same number of "
                "binding sites"<<endl;
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;
        bool state = false;
        auto cyl1_o = get<0>(it1->first)->getCylinder();
        auto bs1_o =  get<1>(it1->first);
        auto cyl2_o = get<0>(it1->second)->getCylinder();
        auto bs2_o =  get<1>(it1->second);
        //
        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(it2->first)->getCylinder();
            auto bs1_s =  get<1>(it2->first);
            auto cyl2_s = get<0>(it2->second)->getCylinder();
            auto bs2_s =  get<1>(it2->second);
            sum = 0;
            if(cyl1_o->getID() == cyl1_s->getID() )
                sum++;
            if(cyl2_o->getID() == cyl2_s->getID() )
                sum++;
            if(bs1_o == bs1_s )
                sum++;
            if(bs2_o == bs2_s )
                sum++;
            if(sum == 4) {
//                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else {
                    cout << "ERROR. Multiple matches for chosen binding site pair in "
                            "STENCILLIST. Check stencillist.." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    std::cout<<"possible bindings size "<<_possibleBindings.size()<<" Total matches "<<
             matches<<endl;
        if(_possibleBindings.size() != matches || _possibleBindings.size() !=
                                                      _possibleBindingsstencil.size()){
        cout<<"Linker. All binding site pairs are not found in Stencil list"<<endl;
            exit(EXIT_FAILURE);
    }
}
#endif
#ifdef CUDAACCL_NL
void LinkerBindingManager::assigncudavars() {
//    if(gpu_rminmax == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_rminmax, 2 * sizeof(double)), "cuda data transfer", " "
                "BindingManager.cu");
        double dist[2];
        dist[0] = _rMin;
        dist[1] = _rMax;
        CUDAcommon::handleerror(cudaMemcpy(gpu_rminmax, dist, 2 *sizeof(double), cudaMemcpyHostToDevice));
//    }
//    if(gpu_numpairs == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)), "cuda data transfer", " "
                "BindingManager.cu");
//    int n[1];
//    n[0] = 0;
//    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
//    }
//    delete dist;
}
void LinkerBindingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_rminmax),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
}

#endif
/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> LinkerBindingManager::chooseBindingSites() {

    assert((_possibleBindings.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
           sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();
    advance(it, randomIndex);
#ifdef DETAILEDOUTPUT
    auto xxx = _compartment->coordinates();
    std::cout<<"Compartment coords "<<xxx[0]<<" "<<xxx[1]<<" "<<xxx[2]<<endl;
#endif
#ifdef DEBUGCONSTANTSEED
    return (*it);
//    return vector<tuple<CCylinder*, short>>{it->first, it->second};
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
//                ""<<c2[1]<<" "<<c2[1]<<" "<<c2[2]<<endl;
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

#ifdef NLORIGINAL
void MotorBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

    if (cc->getType() != _filamentType) return;
#ifdef DEBUGCONSTANTSEED
    struct Orderset
    {
        bool operator()(Cylinder* lhs, Cylinder* rhs) const  {
            return lhs->getID() < rhs->getID();
        }
    };
#endif
    //if we change other managers copy number
    vector<MotorBindingManager *> affectedManagers;

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

        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : nList) {

            Cylinder * c = cc->getCylinder();

            if (cn->getParent() == c->getParent()) continue;
            if (cn->getType() != _filamentType) continue;

            auto ccn = cn->getCCylinder();

            for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {

                if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(),
                             1.0)) {

                    //check distances..
                    auto mp1 = (float) bindingSite /
                               SysParams::Geometry().cylinderNumMon[_filamentType];
                    auto mp2 = (float) *it /
                               SysParams::Geometry().cylinderNumMon[_filamentType];

                    auto x1 = c->getFirstBead()->coordinate;
                    auto x2 = c->getSecondBead()->coordinate;
                    auto x3 = cn->getFirstBead()->coordinate;
                    auto x4 = cn->getSecondBead()->coordinate;

                    auto m1 = midPointCoordinate(x1, x2, mp1);
                    auto m2 = midPointCoordinate(x3, x4, mp2);

                    double dist = twoPointDistance(m1, m2);

                    if (dist > _rMax || dist < _rMin) continue;

                    auto t1 = tuple<CCylinder *, short>(cc, bindingSite);
                    auto t2 = tuple<CCylinder *, short>(ccn, *it);

                    //add in correct order
                    if (c->getID() > cn->getID()) {
#ifdef DEBUGCONSTANTSEED
                        appendpossibleBindings(t1, t2);
//                        _possibleBindings.emplace(t1, t2);
#else
                        _possibleBindings.emplace(t1, t2);
#endif
                    } else {
                        //add in this compartment
                            if (cn->getCompartment() == _compartment) {

#ifdef DEBUGCONSTANTSEED
                                appendpossibleBindings(t1, t2);
//                            _possibleBindings.emplace(t1, t2);
#else
                                _possibleBindings.emplace(t1, t2);
#endif
                            }
                                //add in other
                            else {

                                auto m = (MotorBindingManager *) cn->getCompartment()->
                                        getFilamentBindingManagers()[_mIndex].get();

                                affectedManagers.push_back(m);

#ifdef DEBUGCONSTANTSEED
                                m->appendpossibleBindings(t2, t1);
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
        for (auto m : affectedManagers) {

            int oldNOther = m->_bindingSpecies->getN();
            int newNOther = m->numBindingSites();

            m->updateBindingReaction(oldNOther, newNOther);
        }

        //update this manager
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSites();

        updateBindingReaction(oldN, newN);
    }
#endif

void MotorBindingManager::addPossibleBindings(CCylinder* cc) {
//    auto cylcoord = cc->getCylinder()->coordinate;
//    std::cout<<"M Adding possible to Cyl "<<cc->getCylinder()->getID()<<" w coords "
//            ""<<cylcoord[0]<<" "<<cylcoord[1]<<" "<<cylcoord[2]<<endl;

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
}


void MotorBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
#ifdef DEBUGCONSTANTSEED
    erasepossibleBindings(cc,bindingSite);
    //remove all tuples which have this ccylinder as key
//    auto t = tuple<CCylinder*, short>(cc, bindingSite);
//    _possibleBindings.erase(t);
//
//    //remove all tuples which have this as value
//    for (auto it = _possibleBindings.begin(); it != _possibleBindings.end(); ) {
//
//        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
//            _possibleBindings.erase(it++);
//
//        else ++it;
//    }
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

//        for (auto it = m->_possibleBindings.begin(); it != m->_possibleBindings.end(); ) {
//
//            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
//                m->_possibleBindings.erase(it++);
//
//            else ++it;
//        }
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
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
#ifdef NLSTENCILLIST
        removePossibleBindingsstencil(cc, *bit);
#endif
    }
}

#ifdef NLORIGINAL
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
    for(auto c :_cylinderssorted)
#else
//    std::cout<<_compartment->getCylinders().size()<<endl;
        for(auto c :_compartment->getCylinders())
#endif
    {
//        auto x = _compartment->coordinates();
//        std::cout<<"M updateall Cyl "<<c->getID()<<" "<<c->coordinate[0]<<" "<<
//                 c->coordinate[1]<<" "<<c->coordinate[2]<<" in comp "<<x[0]<<" "
//                ""<<x[1]<<" "<<x[2]<<endl;

        if(c->getType() != _filamentType) continue;

        auto cc = c->getCCylinder();

        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {

            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {

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
                            if(c->getID() > cn->getID()) {
                                //                                        ""<<x[1]<<" "<<x[2]<<endl;
#ifdef DEBUGCONSTANTSEED
                                appendpossibleBindings(t1,t2);
//                                _possibleBindings.emplace(t1, t2);
#else
                                _possibleBindings.emplace(t1, t2);
#endif
//                                std::cout<<"M pb Cyl "<<cc->getCylinder()->getID()<<" bs "
//                                        ""<<*it1<<" Cyl "<<ccn->getCylinder()->getID()<<""
//                                        " bs "<<*it2<<" in comp "<<x[0]<<" "
//                                        ""<<x[1]<<" "<<x[2]<<endl;
                            }
                        }
                    }
                }
            }
        }
    }
//    auto x = _compartment->coordinates();
//    std::cout<<"M possibile "<<_possibleBindings.size()<<" in comp "<<x[0]<<" "
//                     ""<<x[1]<<" "<<x[2]<<endl;
//    std::cout<<"*****"<<endl;
//    for(auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++){
//        auto bspair = *it;
//        auto cyl1 = get<0>(bspair[0])->getCylinder()->getID();
//        auto cyl2 = get<0>(bspair[1])->getCylinder()->getID();
//        std::cout<<"Cyl "<<cyl1<<"bs "<<get<1>(bspair[0])<<" Cyl "<<cyl2<<" bs "
//                ""<<get<1>(bspair[1])<<endl;
//    }
//    for(auto it = _possibleBindings.begin(); it != _possibleBindings.end(); it++){
//        auto cyl1 = get<0>(it->first)->getCylinder();
//        auto cyl2 = get<0>(it->second)->getCylinder();
//        std::cout<<"Cyl "<<cyl1->getID()<<"bs "<<get<1>(it->first)<<" Cyl "<<cyl2->getID
//                ()<<" bs "<<get<1>(it->second)<<endl;
//    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
//std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<_compartment->coordinates()[1]<<""
//        " "<<_compartment->coordinates()[2]<<" "<<_possibleBindings.size()<<" "<<oldN<<" "<<newN<<endl;
    updateBindingReaction(oldN, newN);
}
#endif

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
#ifdef NLSTENCILLIST
void MotorBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {

        //loop through neighbors
        //now re add valid based on CCNL
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

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
                        _possibleBindingsstencil.emplace(t1,t2);
                    }

                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {

                            _possibleBindingsstencil.emplace(t2,t1);
                        }
                            //add in other
                        else {

                            auto m = (MotorBindingManager*)cn->getCompartment()->
                                    getFilamentBindingManagers()[_mIndex].get();

                            affectedManagers.push_back(m);

                            m->_possibleBindingsstencil.emplace(t2,t1);
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
void MotorBindingManager::updateAllPossibleBindingsstencil() {

    _possibleBindingsstencil.clear();

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
                for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil
                        (cc->getCylinder())) {

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
                                _possibleBindingsstencil.emplace(t1,t2);
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
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)

        removePossibleBindingsstencil(cc, *bit);
}
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencil.erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencil.begin(); it != _possibleBindingsstencil.end();
            ) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindingsstencil.erase(it++);

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

void MotorBindingManager::crosscheck(){
    //cout<<"Branching NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
    //    <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size()) {
        cout << "Motor.. The two methods compared do not yield the same number of "
                "binding sites" << endl;
        exit(EXIT_FAILURE);
    }
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;
        bool state = false;
        auto cyl1_o = get<0>(it1->first)->getCylinder();
        auto bs1_o =  get<1>(it1->first);
        auto cyl2_o = get<0>(it1->second)->getCylinder();
        auto bs2_o =  get<1>(it1->second);
        //
        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(it2->first)->getCylinder();
            auto bs1_s =  get<1>(it2->first);
            auto cyl2_s = get<0>(it2->second)->getCylinder();
            auto bs2_s =  get<1>(it2->second);
            sum = 0;
            if(cyl1_o->getID() == cyl1_s->getID() )
                sum++;
            if(cyl2_o->getID() == cyl2_s->getID() )
                sum++;
            if(bs1_o == bs1_s )
                sum++;
            if(bs2_o == bs2_s )
                sum++;
            if(sum == 4) {
//                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else
                    cout<<"ERROR. Multiple matches for chosen binding site pair in "
                            "STENCILLIST. Check stencillist.."<<endl;
            }
        }
    }
    std::cout<<"possible bindings size "<<_possibleBindings.size()<<" Total matches "<<
             matches<<endl;
    if(_possibleBindings.size() != matches || _possibleBindings.size() !=
                                              _possibleBindingsstencil.size()){
        cout<<"Motor. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
#endif
#ifdef CUDAACCL_NL
void MotorBindingManager::assigncudavars() {

//    if(gpu_numpairs == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)), "cuda data transfer", " "
                "BindingManager.cu");
//    int n[1];
//    n[0] = 0;
//    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
//    }

//    if(gpu_rminmax == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_rminmax, 2 * sizeof(double)), "cuda data transfer", " "
                "BindingManager.cu");
        double dist[2];
        dist[0] = _rMin;
        dist[1] = _rMax;
        CUDAcommon::handleerror(cudaMemcpy(gpu_rminmax, dist, 2 * sizeof(double), cudaMemcpyHostToDevice));
//    }

//    delete dist;
}
void MotorBindingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_rminmax),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
}
#endif


/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> MotorBindingManager::chooseBindingSites() {

    assert((_possibleBindings.size() != 0)
           && "Major bug: Motor binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
    auto it = _possibleBindings.begin();

    advance(it, randomIndex);
//    auto xxx = _compartment->coordinates();
//    std::cout<<"Compartment coords "<<xxx[0]<<" "<<xxx[1]<<" "<<xxx[2]<<endl;

#ifdef DEBUGCONSTANTSEED
    return (*it);
//    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#else
    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}
#ifdef DEBUGCONSTANTSEED
    void MotorBindingManager::erasepossibleBindings(CCylinder* cc, short bindingSite) {
        tuple<CCylinder*, short> bindingtoremove = make_tuple(cc, bindingSite);
//        std::cout<<"erasing cyl "<<cc->getCylinder()->getID()<<" bs "<<bindingSite<<endl;
        int counter = 0;
        for (auto p = _possibleBindings.begin(); p != _possibleBindings.end(); p++) {
            auto binding1 = (*p)[0];
            auto binding2 = (*p)[1];
            if (bindingtoremove == binding1 || bindingtoremove == binding2)
            {
//                auto cyl1 = get<0>(binding1)->getCylinder()->getID();
//                auto cyl2 = get<0>(binding2)->getCylinder()->getID();
//                auto bs1 = get<1>(binding1);
//                auto bs2 = get<1>(binding2);
//                auto x = _compartment->coordinates();
//                std::cout<<"Removing pair Cyl "<<cyl1<<" bs "<<bs1<<" Cyl "<<cyl2<<" bs "
//                        ""<<bs2<<"  ID "<<counter<<" in compartment "<<x[0]<<" "<<x[1]<<" "
//                                 ""<<x[2]<<endl;
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
//            auto c1 = (get<0>(t1))->getCylinder()->coordinate;
//            auto c2 = (get<0>(t2))->getCylinder()->coordinate;
//            auto c1ID = (get<0>(t1))->getCylinder()->getID();
//            auto c2ID = (get<0>(t2))->getCylinder()->getID();
//            auto bs1 = get<1>(t1);
//            auto bs2 = get<1>(t2);
//            auto x = _compartment->coordinates();
//            std::cout<<"Added Cyl "<<c1ID<<" bs "<<bs1<<" Cyl "<<c2ID<<" bs "<<bs2<<" "
//                    "to position "<<_possibleBindings.size()<<" in compartment "<<x[0]<<" "<<x[1]<<" "
//                             ""<<x[2]<<"with coords "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "
//                             ""<<c2[1]<<" "<<c2[1]<<" "<<c2[2]<<endl;
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
//        std::cout<<"oldN "<<oldN<<" newN "<<newN<<endl;
        updateBindingReaction(oldN,newN);
    }

SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CylinderCylinderNL*> LinkerBindingManager::_neighborLists;
vector<CylinderCylinderNL*> MotorBindingManager::_neighborLists;

