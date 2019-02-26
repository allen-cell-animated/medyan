
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
#include "HybridBindingSearchManager.h"
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
    }
}

void BranchingManager::updateAllPossibleBindings() {


    //clear all
    _possibleBindings.clear();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    //int offset = SysParams::Mechanics().bsoffsetvec.at(_filamentType);

    for(auto &c : _compartment->getCylinders()) {

        if(c->getType() != _filamentType) continue;

        auto cc = c->getCCylinder();
        int j = -1;
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
            j++;
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
            if (areEqual(boundstate[0][offset + SysParams::Chemistry()
                                       .bindingSites[_filamentType]
                                       .size()*c->_dcIndex + j], 1.0) && inZone) {
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
    /*std::cout<<"Branching consistency "<<isConsistent()<<endl;*/
}
#endif
bool BranchingManager::isConsistent() {
#ifdef NLORIGINAL
    auto bindinglist = _possibleBindings;
#elif defined(NLSTENCILLIST)
    auto bindinglist = _possibleBindingsstencil;
#endif
    for (auto it = bindinglist.begin(); it != bindinglist.end(); it++) {

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
void BranchingManager::addPossibleBindingsstencil(CCylinder* cc) {
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
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
            else if(_nucleationZone == NucleationZoneType::RIGHTBOUNDARY){
                double dis = _subSystem->getBoundary()->getboundaryelementcoord(1) -
                        coord[0];
                if(dis > 0 && dis < _nucleationDistance)
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
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
}
void BranchingManager::updateAllPossibleBindingsstencil() {
    //clear all
    _possibleBindingsstencil.clear();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    //int offset = SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    for(auto &c : _compartment->getCylinders()) {

        if(c->getType() != _filamentType) continue;

        auto cc = c->getCCylinder();
        int j = -1;
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
            j++;
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
                    else if(_nucleationZone == NucleationZoneType::RIGHTBOUNDARY){
                        double dis = _subSystem->getBoundary()->getboundaryelementcoord(1) -
                                     coord[0];
                        if(dis > 0 && dis < _nucleationDistance)
                            inZone = true;
                        else
                            inZone = false;
                    }

                    else inZone = true;
                }
                else
                inZone = false;
            }
            if (areEqual(boundstate[0][maxnbs * c->_dcIndex + j], 1.0) && inZone) {
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
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
    /*std::cout<<"Branching consistency "<<isConsistent()<<endl;*/
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
    int newN = numBindingSitesstencil();

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
    std::cout<<"Branching possible bindings size "<<_possibleBindings.size()<<" Total "
            "matches "<<
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
    _rMinsq =_rMin * _rMin;
    _rMaxsq = _rMax * _rMax;

    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[ML_RXN_INDEX]->getSpecies().getName();

    _bindingSpecies = _compartment->findSpeciesByName(name);
    _rMaxsq = rMax*rMax;
    _rMinsq = rMin*rMin;
    minparamcyl2 = (float)*(SysParams::Chemistry().bindingSites[_filamentType].begin())/
                   SysParams::Geometry().cylinderNumMon[_filamentType];
    maxparamcyl2 = (float)(SysParams::Chemistry().bindingSites[_filamentType].back())/
                   SysParams::Geometry().cylinderNumMon[_filamentType];
    for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
        it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        bindingsites.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[_filamentType]);
    }
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

                    double distsq = twoPointDistancesquared(m1, m2);

                    if(distsq > _rMaxsq || distsq < _rMinsq) continue;

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
    }
}

void LinkerBindingManager::updateAllPossibleBindings() {

    _possibleBindings.clear();
    double min1,min2,max1,max2;
    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;
    double timetaken = 0.0;
    double time16 = 0.0;
    double minparamcyl2 = (float)*(SysParams::Chemistry().bindingSites[_filamentType].begin())/
    SysParams::Geometry().cylinderNumMon[_filamentType];
    double maxparamcyl2 = (float)(SysParams::Chemistry().bindingSites[_filamentType].back())/
    SysParams::Geometry().cylinderNumMon[_filamentType];
    double sqdisttermswithjustalpha;
    bool status1 = true;
    bool status2 = true;
    vector<double> maxvec;
    vector<double> minvec;
    int accepts = 0;
    int total = 0;
    int rejects16 = 0;
    int rejectsnavail =0;
    mints = chrono::high_resolution_clock::now();
    vector<double> bindingsites;
    vector<double> cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;

    for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
        it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        bindingsites.push_back((float)*it1 / SysParams::Geometry()
                               .cylinderNumMon[_filamentType]);
    }

    //    minte = chrono::high_resolution_clock::now();
    //    chrono::duration<double> elapsed_vec(minte - mints);
    //    std::cout<<"Vectorize time "<<elapsed_vec.count()<<endl;

    accepts =0;
    total = 0;
    time16 = 0.0;
    timetaken = 0.0;
    _possibleBindings.clear();
    mints = chrono::high_resolution_clock::now();
    //int offset = SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    for(auto c : _compartment->getCylinders()) {
        if (c->getType() != _filamentType) continue;

        auto x1 = c->getFirstBead()->coordinate;
        auto x2 = c->getSecondBead()->coordinate;
        auto cc = c->getCCylinder();
        vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
            if(c->getID() < cn->getID()) continue;
            auto ccn = cn->getCCylinder();
            auto x3 = cn->getFirstBead()->coordinate;
            auto x4 = cn->getSecondBead()->coordinate;

            vector<double> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            vector<double> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            double maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);

            double mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
                                                                                    x3,x4)));
            mindistsq = mindistsq * mindistsq;
            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;

            double X1X3squared = sqmagnitude(X1X3);
            double X1X2squared = cylsqmagnitudevector.at(c->_dcIndex);
            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            double X3X4squared = cylsqmagnitudevector.at(cn->_dcIndex);
            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            mins2 = chrono::high_resolution_clock::now();
            int i = -1;
            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                i++;
                //now re add valid binding sites
                if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                           .bindingSites[_filamentType].size()
                                           *c->_dcIndex + i], 1.0)) {
                    auto mp1 = bindingsites.at(i);
                    double A = X3X4squared;
                    double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 *
                    X1X3dotX1X2;
                    double C1 = C - _rMinsq;
                    double C2 = C - _rMaxsq;
                    double b2m4ac1 = B*B - 4*A*C1;
                    double b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) continue;
                    maxvec.clear();
                    minvec.clear();
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minvec.push_back(min1);
                            minvec.push_back(min2);
                        }
                        else{
                            minvec.push_back(min2);
                            minvec.push_back(min1);
                        }
                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxvec.push_back(max1);
                            maxvec.push_back(max2);
                        }
                        else{
                            maxvec.push_back(max2);
                            maxvec.push_back(max1);
                        }
                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
                            continue;
                        }
                    }
                    int j =-1;
                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                        it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
                        j++;
                        bool check2 = true;
                        if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                                   .bindingSites[_filamentType]
                                                   .size()*cn->_dcIndex + j], 1.0)) {
                            total++;
                            //check distances..
                            auto mp2 = bindingsites.at(j);

                            if(!status2) {
                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1)) {
                                    //                                    check2 = false;
                                    continue;
                                }
                            }
                            if(!status1){
                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1)) {
                                    //                                    check2 = false;
                                    continue;
                                }
                            }

                            accepts++;
                            if(check2) {
                                mins = chrono::high_resolution_clock::now();
                                auto t1 = tuple<CCylinder *, short>(cc, *it1);
                                auto t2 = tuple<CCylinder *, short>(ccn, *it2);

                                //add in correct order
                                if (c->getID() > cn->getID()) {
                                    _possibleBindings.emplace(t1, t2);
                                }
                            }
                            mine= chrono::high_resolution_clock::now();
                            chrono::duration<double> elapsed_emplace(mine - mins);
                            timetaken += elapsed_emplace.count();
                        }
                    }
                }
            }
            mine2= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_run16(mine2 - mins2);
            time16 += elapsed_run16.count();
        }
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    updateBindingReaction(oldN, newN);
}

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
#endif
bool LinkerBindingManager::isConsistent() {

#ifdef NLORIGINAL
    auto bindinglist = _possibleBindings;
#elif defined(NLSTENCILLIST)
    auto bindinglist = _possibleBindingsstencil;
#endif
    for (auto it = bindinglist.begin(); it != bindinglist.end(); it++) {
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
#ifdef NLSTENCILLIST
void LinkerBindingManager::addPossibleBindingsstencil(CCylinder* cc) {
/*    std::cout<<"Adding possible bindings of cylinder with ID "<<cc->getCylinder()
            ->getID()<<" with cindex "<<cc->getCylinder()->_dcIndex<<endl;*/
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
void LinkerBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#ifdef HYBRID_NLSTENCILLIST
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->addPossibleBindingsstencil(_idvec,cc,bindingSite);
#else
    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                                                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0)) {

        //loop through neighbors
        //now re add valid based on CCNL
        vector<Cylinder*> Neighbors;
//#ifdef HYBRID_NLSTENCILLIST
//        Neighbors = cc->getCompartment()->getHybridBindingSearchManager()->getHNeighbors
//                (cc->getCylinder(),HNLID);
//#else
        Neighbors = _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder());
//#endif
        for (auto cn : Neighbors) {

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

                    double distsq = twoPointDistancesquared(m1, m2);

                    if(distsq > _rMaxsq || distsq < _rMinsq) continue;

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
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }

    //update this manager
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
#endif
}
void LinkerBindingManager::updateAllPossibleBindingsstencil() {
    _possibleBindingsstencil.clear();
    //int offset = SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    double min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    double minveca[2];
    double maxveca[2];
    int accepts = 0;
    int total = 0;
    int Ncyl = Cylinder::getCylinders().size();
    int nbs = SysParams::Chemistry().bindingSites[_filamentType].size();
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    double* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    double* coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    CCylinder** ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    auto bindingsitevec =SysParams::Chemistry().bindingSites[_filamentType];
    int Ncylincmp =  _compartment->getCylinders().size();
    int cindexvec[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int>ncindex; //helper vector
    long id = 0;
    for(auto c : _compartment->getCylinders()){
        cindexvec[id] = c->_dcIndex;
//        totalneighbors += _neighborLists[_nlIndex]->getNeighborsstencil(c).size();
        id++;
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(c)) {
            if(c->getID() > cn->getID())
                ncindex.push_back(cn->_dcIndex);
//            else
//                counter1++;
        }
        ncindices.push_back(ncindex);
        ncindex.clear();
    }
    //Check 1
/*    for(int idx = 0; idx < Ncyl; idx++) {
        auto a1 = cylindervec[idx].cindex;
        auto a2 = ccylvec[a1]->getCylinder()->_dcIndex;
        if(a1 != a2)
            std::cout <<"Cyl ID mismatch M "<<cylindervec[idx].cindex << " " <<
                      ccylvec[a1]->getCylinder()->_dcIndex << endl;
    }*/

    for(int i=0;i<Ncylincmp;i++){
        int cindex = cindexvec[i];
        cylinder c = cylindervec[cindex];
        if(c.type != _filamentType) {
//            counter2++;
            continue;}
        double x1[3],x2[3];
        memcpy(x1, &coord[3*c.bindices[0]], 3 * sizeof(double));
        memcpy(x2, &coord[3*c.bindices[1]], 3 * sizeof(double));
        double X1X2[3] ={x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        //Check 2
/*        for(auto cndummy:ncindices[i]){
            auto A =cylindervec[cndummy];
            if(ccylvec[A.cindex]->getCylinder()->getID() != A.ID)
            {
                std::cout<<"Mismatch in neighbors L of Cyl index "<<i<<" "
                        ""<<ccylvec[cndummy]->getCylinder()->getID()<<" "
                        ""<<cylindervec[cndummy].ID<<endl;
            }
        }*/

        int* cnindices = ncindices[i].data();
        for(int arraycount = 0; arraycount < ncindices[i].size();arraycount++){
            int cnindex = cnindices[arraycount];
//            std::cout<<*cnindex<<endl;
            cylinder cn = cylindervec[cnindex];
//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
// not contain ncs that will fail this cndn.
            if(c.filamentID == cn.filamentID){
//                counter3++;
                continue;}
            if(c.type != cn.type){
//                counter4++;
                continue;}

            double x3[3], x4[3];
            memcpy(x3, &coord[3*cn.bindices[0]], 3 * sizeof(double));
            memcpy(x4, &coord[3*cn.bindices[1]], 3 * sizeof(double));
            double X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            double X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            double X1X3squared = sqmagnitude(X1X3);
            double X1X2squared = cylsqmagnitudevector[c.cindex];
            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            double X3X4squared = cylsqmagnitudevector[cn.cindex];
            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            for(int pos1 =0; pos1<nbs;pos1++) {
                //now re add valid binding sites
                if (areEqual(boundstate[1][maxnbs *c.cindex + pos1], 1.0)) {
                    auto mp1 = bindingsites.at(pos1);
                    double A = X3X4squared;
                    double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 * X1X3dotX1X2;
                    double C1 = C - _rMinsq;
                    double C2 = C - _rMaxsq;
                    double b2m4ac1 = B*B - 4*A*C1;
                    double b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) {
//                    counter6++;
                        continue;}
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minveca[0] = (min1);
                            minveca[1] = (min2);
                        }
                        else{
                            minveca[0] = (min2);
                            minveca[1] = (min1);
                        }
                        if(minveca[0]< minparamcyl2 && minveca[1] > maxparamcyl2) {
//                        counter7++;
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxveca[0] = (max1);
                            maxveca[1] = (max2);
                        }
                        else{
                            maxveca[0] = (max2);
                            maxveca[1] = (max1);
                        }
                        if(maxveca[0] > maxparamcyl2 || maxveca[1] <minparamcyl2){
//                        counter8++;
                            continue;
                        }
                    }
                    for(int pos2 = 0; pos2<nbs;pos2++){
                        if (areEqual(boundstate[1][maxnbs *cn.cindex + pos2], 1.0)) {
                            total++;
                            //check distances..
                            auto mp2 = bindingsites.at(pos2);
                            if(!status2) {
                                if (mp2 < maxveca[0] || mp2 > maxveca[1]) {
                                    {
//                                    counter9++;
                                        continue;}
                                }
                            }
                            if(!status1){
                                if (mp2 > minveca[0] && mp2 < minveca[1]) {
                                    {
//                                    counter10++;
                                        continue;}
                                }
                            }
                            accepts++;
                            auto it1 = SysParams::Chemistry().bindingSites[_filamentType][pos1];
                            auto it2 = SysParams::Chemistry().bindingSites[_filamentType][pos2];
                            auto xx1 = ccylvec[cindex]->getCylinder();
                            auto xx2 = ccylvec[cnindex]->getCylinder();
//                            std::cout<<xx1->getID()<<" "<<c.ID<<" "<<xx2->getID()<<" "
//                                    ""<<cn.ID<<endl;
                            /*if(xx1->getFirstBead()->_dbIndex != c.bindices[0] ||
                               xx1->getSecondBead()->_dbIndex != c.bindices[1])
                                std::cout<<"DB1 "<<xx1->getFirstBead()->_dbIndex<<" "
                                        ""<<xx1->getSecondBead()->_dbIndex<<" "<<c
                                                 .bindices[0]<<" "<<c.bindices[1]<<endl;
                            if(xx2->getFirstBead()->_dbIndex != cn.bindices[0] ||
                               xx2->getSecondBead()->_dbIndex != cn.bindices[1])
                                std::cout<<"DB2 "<<xx2->getFirstBead()->_dbIndex<<" "
                                        ""<<xx2->getSecondBead()->_dbIndex<<" "<<cn
                                                 .bindices[0]<<" "<<cn.bindices[1]<<endl;*/
                            auto t1 = tuple<CCylinder *, short>(ccylvec[cindex], it1);
                            auto t2 = tuple<CCylinder *, short>(ccylvec[cnindex], it2);
                            //add in correct order
                            _possibleBindingsstencil.emplace(t1, t2);
                        }
                    }
                }
            }
        }
    }
//    std::cout<<"Tuple size "<<_possibleBindingsstencil.size()<<endl;
    /*_possibleBindingsstencil.clear();
    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;
    double timetaken = 0.0;
    double time16 = 0.0;
    double sqdisttermswithjustalpha;
    vector<double> maxvec;
    vector<double> minvec;
    int rejects16 = 0;
    int rejectsnavail =0;
    mints = chrono::high_resolution_clock::now();
    vector<double> bindingsites;

    for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
        it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        bindingsites.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[_filamentType]);
    }

    minte = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_vec(minte - mints);
    std::cout<<"Vectorize time "<<elapsed_vec.count()<<endl;

    accepts =0;
    total = 0;
    time16 = 0.0;
    timetaken = 0.0;
    _possibleBindingsstencil.clear();
    mints = chrono::high_resolution_clock::now();
    for(auto c : _compartment->getCylinders()) {
        if (c->getType() != _filamentType) continue;
        auto x1 = c->getFirstBead()->coordinate;
        auto x2 = c->getSecondBead()->coordinate;
        auto cc = c->getCCylinder();
        vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
            if(c->getID() < cn->getID()) continue;
            auto ccn = cn->getCCylinder();
            auto x3 = cn->getFirstBead()->coordinate;
            auto x4 = cn->getSecondBead()->coordinate;

            vector<double> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            vector<double> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            double maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);

            double mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
                                                                                    x3,x4)));
            mindistsq = mindistsq * mindistsq;
            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;

            double X1X3squared = sqmagnitude(X1X3);
            double X1X2squared = cylsqmagnitudevector.at(c->_dcIndex);
            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            double X3X4squared = cylsqmagnitudevector.at(cn->_dcIndex);
            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            mins2 = chrono::high_resolution_clock::now();
            int i = -1;
            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                i++;
                //now re add valid binding sites
                if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                                            .bindingSites[_filamentType].size()
                                                    *c->_dcIndex + i], 1.0)) {
                    auto mp1 = bindingsites.at(i);
                    double A = X3X4squared;
                    double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 *
                                                                       X1X3dotX1X2;
                    double C1 = C - _rMinsq;
                    double C2 = C - _rMaxsq;
                    double b2m4ac1 = B*B - 4*A*C1;
                    double b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) continue;
                    maxvec.clear();
                    minvec.clear();
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minvec.push_back(min1);
                            minvec.push_back(min2);
                        }
                        else{
                            minvec.push_back(min2);
                            minvec.push_back(min1);
                        }
                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxvec.push_back(max1);
                            maxvec.push_back(max2);
                        }
                        else{
                            maxvec.push_back(max2);
                            maxvec.push_back(max1);
                        }
                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
                            continue;
                        }
                    }
                    int j =-1;
                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                        it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
                        j++;
                        bool check2 = true;
                        if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                                                    .bindingSites[_filamentType]
                                                                    .size()*cn->_dcIndex + j], 1.0)) {
                            total++;
                            //check distances..
                            auto mp2 = bindingsites.at(j);

                            if(!status2) {
                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1)) {
//                                    check2 = false;
                                    continue;
                                }
                            }
                            if(!status1){
                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1)) {
//                                    check2 = false;
                                    continue;
                                }
                            }

                            accepts++;
                            if(check2) {
                                mins = chrono::high_resolution_clock::now();
                                auto t1 = tuple<CCylinder *, short>(cc, *it1);
                                auto t2 = tuple<CCylinder *, short>(ccn, *it2);

                                //add in correct order
                                if (c->getID() > cn->getID()) {
                                    _possibleBindingsstencil.emplace(t1, t2);
                                }
                            }
                            mine= chrono::high_resolution_clock::now();
                            chrono::duration<double> elapsed_emplace(mine - mins);
                            timetaken += elapsed_emplace.count();
                        }
                    }
                }
            }
            mine2= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_run16(mine2 - mins2);
            time16 += elapsed_run16.count();
        }
    }
    minte = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_total3(minte - mints);*/
//    std::cout<<"Overall time "<<elapsed_total3.count()<<endl;
//    std::cout<<"Total "<<total<<" accepts "<<accepts<<endl;
//    std::cout<<"16 loop time taken "<<time16<<endl;
//    std::cout<<"Tuple Emplace time taken "<<timetaken<<endl;
//    std::cout<<"Tuple size "<<_possibleBindingsstencil.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();
    updateBindingReaction(oldN, newN);
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        removePossibleBindingsstencil(cc, *bit);
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#ifdef HYBRID_NLSTENCILLIST
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->removePossibleBindingsstencil(_idvec, cc, bindingSite);
/*    for(auto C:SubSystem::getstaticgrid()->getCompartments()){
        C->getHybridBindingSearchManager()->checkoccupancySIMD(_idvec);
    }*/
#else

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
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);

    //remove all neighbors which have this binding site pair
/*#ifdef HYBRID_NLSTENCILLIST
    auto C = cc->getCompartment();
    for (auto cn : C->getHybridBindingSearchManager()->getHNeighbors(cc->getCylinder(), HNLID)){
        if(cn->getType() != _filamentType) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (LinkerBindingManager*)cn->getCompartment()->
                    getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
#else*/
    for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

        if(cn->getType() != _filamentType) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (LinkerBindingManager*)cn->getCompartment()->
            getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
            affectedManagers.push_back(m);
        }
    }
//#endif

    //remove, update affected
    for(auto m : affectedManagers) {

        for (auto it = m->_possibleBindingsstencil.begin(); it !=
             m->_possibleBindingsstencil.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            m->_possibleBindingsstencil.erase(it++);

            else ++it;
        }

        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }
#endif
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
    std::cout<<"Linker possible bindings size "<<_possibleBindings.size()<<" Total matches"
            " "<<
             matches<<endl;
        if(_possibleBindings.size() != matches || _possibleBindings.size() !=
                                                      _possibleBindingsstencil.size()){
        cout<<"Linker. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
vector<tuple<CCylinder*, short>> LinkerBindingManager::chooseBindingSitesstencil() {
#ifdef HYBRID_NLSTENCILLIST
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->chooseBindingSitesstencil(_idvec);
#else
    assert((_possibleBindingsstencil.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindingsstencil.size() - 1);
    auto it = _possibleBindingsstencil.begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
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
//MOTOR
MotorBindingManager::MotorBindingManager(ReactionBase* reaction,
                                         Compartment* compartment,
                                         short boundInt, string boundName,
                                         short filamentType,
                                         float rMax, float rMin)

: FilamentBindingManager(reaction, compartment, boundInt, boundName, filamentType),
_rMin(rMin), _rMax(rMax) {

    _rMinsq =_rMin * _rMin;
    _rMaxsq = _rMax * _rMax;

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
    _rMaxsq = rMax*rMax;
    _rMinsq = rMin*rMin;
    minparamcyl2 = (float)*(SysParams::Chemistry().bindingSites[_filamentType].begin())/
                          SysParams::Geometry().cylinderNumMon[_filamentType];
    maxparamcyl2 = (float)(SysParams::Chemistry().bindingSites[_filamentType].back())/
                          SysParams::Geometry().cylinderNumMon[_filamentType];
    for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
        it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        bindingsites.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[_filamentType]);
    }
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

                    double distsq = twoPointDistancesquared(m1, m2);

                    if(distsq > _rMaxsq || distsq < _rMinsq) continue;

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
    }
}

void MotorBindingManager::updateAllPossibleBindings() {

    _possibleBindings.clear();
    double min1,min2,max1,max2;
    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;
    double timetaken = 0.0;
    double time16 = 0.0;
    double minparamcyl2 = (float)*(SysParams::Chemistry().bindingSites[_filamentType].begin())/
    SysParams::Geometry().cylinderNumMon[_filamentType];
    double maxparamcyl2 = (float)(SysParams::Chemistry().bindingSites[_filamentType].back())/
    SysParams::Geometry().cylinderNumMon[_filamentType];
    double sqdisttermswithjustalpha;
    bool status1 = true;
    bool status2 = true;
    vector<double> maxvec;
    vector<double> minvec;
    int accepts = 0;
    int total = 0;
    int rejects16 = 0;
    int rejectsnavail =0;
    mints = chrono::high_resolution_clock::now();
    vector<double> bindingsites;
    vector<double> cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;

    for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
        it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
        bindingsites.push_back((float)*it1 / SysParams::Geometry()
                               .cylinderNumMon[_filamentType]);
    }

    //    minte = chrono::high_resolution_clock::now();
    //    chrono::duration<double> elapsed_vec(minte - mints);
    //    std::cout<<"Vectorize time "<<elapsed_vec.count()<<endl;

    accepts =0;
    total = 0;
    time16 = 0.0;
    timetaken = 0.0;
    _possibleBindings.clear();
    mints = chrono::high_resolution_clock::now();
    //int offset = SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    for(auto c : _compartment->getCylinders()) {
        if (c->getType() != _filamentType) continue;

        auto x1 = c->getFirstBead()->coordinate;
        auto x2 = c->getSecondBead()->coordinate;
        auto cc = c->getCCylinder();
        vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != _filamentType) continue;
            if(c->getID() < cn->getID()) continue;
            auto ccn = cn->getCCylinder();
            auto x3 = cn->getFirstBead()->coordinate;
            auto x4 = cn->getSecondBead()->coordinate;

            vector<double> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            vector<double> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            double maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);

            double mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
                                                                                    x3,x4)));
            mindistsq = mindistsq * mindistsq;
            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;

            double X1X3squared = sqmagnitude(X1X3);
            double X1X2squared = cylsqmagnitudevector.at(c->_dcIndex);
            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            double X3X4squared = cylsqmagnitudevector.at(cn->_dcIndex);
            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            mins2 = chrono::high_resolution_clock::now();
            int i = -1;
            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                i++;
                //now re add valid binding sites
                if (areEqual(boundstate[2][offset + SysParams::Chemistry()
                                           .bindingSites[_filamentType].size()
                                           *c->_dcIndex + i], 1.0)) {
                    auto mp1 = bindingsites.at(i);
                    double A = X3X4squared;
                    double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 *
                    X1X3dotX1X2;
                    double C1 = C - _rMinsq;
                    double C2 = C - _rMaxsq;
                    double b2m4ac1 = B*B - 4*A*C1;
                    double b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) continue;
                    maxvec.clear();
                    minvec.clear();
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minvec.push_back(min1);
                            minvec.push_back(min2);
                        }
                        else{
                            minvec.push_back(min2);
                            minvec.push_back(min1);
                        }
                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxvec.push_back(max1);
                            maxvec.push_back(max2);
                        }
                        else{
                            maxvec.push_back(max2);
                            maxvec.push_back(max1);
                        }
                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
                            continue;
                        }
                    }
                    int j =-1;
                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                        it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
                        j++;
                        bool check2 = true;
                        if (areEqual(boundstate[2][offset + SysParams::Chemistry()
                                                   .bindingSites[_filamentType]
                                                   .size()*cn->_dcIndex + j], 1.0)) {
                            total++;
                            //check distances..
                            auto mp2 = bindingsites.at(j);

                            if(!status2) {
                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1)) {
                                    //                                    check2 = false;
                                    continue;
                                }
                            }
                            if(!status1){
                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1)) {
                                    //                                    check2 = false;
                                    continue;
                                }
                            }

                            accepts++;
                            if(check2) {
                                mins = chrono::high_resolution_clock::now();
                                auto t1 = tuple<CCylinder *, short>(cc, *it1);
                                auto t2 = tuple<CCylinder *, short>(ccn, *it2);

                                //add in correct order
                                if (c->getID() > cn->getID()) {
                                    _possibleBindings.emplace(t1, t2);
                                }
                            }
                            mine= chrono::high_resolution_clock::now();
                            chrono::duration<double> elapsed_emplace(mine - mins);
                            timetaken += elapsed_emplace.count();
                        }
                    }
                }
            }
            mine2= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_run16(mine2 - mins2);
            time16 += elapsed_run16.count();
        }
    }
    //    minte = chrono::high_resolution_clock::now();
    //    chrono::duration<double> elapsed_total3(minte - mints);
    //    std::cout<<"Overall time "<<elapsed_total3.count()<<endl;
    //    std::cout<<"Total "<<total<<" accepts "<<accepts<<endl;
    //    std::cout<<"16 loop time taken "<<time16<<endl;
    //    std::cout<<"Tuple Emplace time taken "<<timetaken<<endl;
    //    std::cout<<"Tuple size "<<_possibleBindings.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    updateBindingReaction(oldN, newN);
    /*std::cout<<"Motor consistency "<<isConsistent()<<endl;*/
}

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
        }
#else
    _possibleBindings.emplace(t1,t2);
#endif
    double newN=numBindingSites();
    updateBindingReaction(oldN,newN);
}
#endif

bool MotorBindingManager::isConsistent() {

#ifdef NLORIGINAL
    auto bindinglist = _possibleBindings;
#elif defined(NLSTENCILLIST)
    auto bindinglist = _possibleBindingsstencil;
#endif
    for (auto it = bindinglist.begin(); it != bindinglist.end(); it++) {

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
void MotorBindingManager::addPossibleBindingsstencil(CCylinder* cc) {
/*    std::cout<<"Adding possible bindings of cylinder with ID "<<cc->getCylinder()
            ->getID()<<" with cindex "<<cc->getCylinder()->_dcIndex<<endl;*/
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
void MotorBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#ifdef HYBRID_NLSTENCILLIST
//    cout<<"Adding "<<cc->getCylinder()->getID()<<" "<<bindingSite<<endl;
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->addPossibleBindingsstencil(_idvec,cc,bindingSite);
//    HManager->checkoccupancySIMD(_idvec);
#else
    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //add valid binding sites
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                                                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {

        //loop through neighbors
        //now re add valid based on CCNL
        vector<Cylinder*> Neighbors;

        Neighbors = _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder());
        for (auto cn : Neighbors) {

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

                    double distsq = twoPointDistancesquared(m1, m2);

                    if(distsq > _rMaxsq || distsq < _rMinsq) continue;

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
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }

    //update this manager
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
#endif
}
void MotorBindingManager::updateAllPossibleBindingsstencil() {

    _possibleBindingsstencil.clear();
    int offset = 0;
            //SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    double min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    double minveca[2];
    double maxveca[2];
    int accepts = 0;
    int total = 0;
    int Ncyl = Cylinder::getCylinders().size();
    int nbs = SysParams::Chemistry().bindingSites[_filamentType].size();
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    double* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    double* coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    CCylinder** ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int counter1 = 0;
    int totalneighbors = 0;
/*    vector<CCylinder*> ccylindervector;
    for(auto cyl:Cylinder::getCylinders())
        ccylindervector.push_back(cyl->getCCylinder());*/

    auto bindingsitevec =SysParams::Chemistry().bindingSites[_filamentType];
    int Ncylincmp =  _compartment->getCylinders().size();
    int cindexvec[Ncylincmp];
    vector<vector<int>> ncindices;
    vector<int>ncindex;
    long id = 0;
    for(auto c : _compartment->getCylinders()){
        cindexvec[id] = c->_dcIndex;
        totalneighbors += _neighborLists[_nlIndex]->getNeighborsstencil(c).size();
        id++;
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(c)) {
            if(c->getID() > cn->getID())
                ncindex.push_back(cn->_dcIndex);
            else
                counter1++;
        }
        ncindices.push_back(ncindex);
        ncindex.clear();
    }
    //Check 1
/*    for(int idx = 0; idx < Ncyl; idx++) {
        auto a1 = cylindervec[idx].cindex;
        auto a2 = ccylvec[a1]->getCylinder()->_dcIndex;
        if(a1 != a2)
            std::cout <<"Cyl ID mismatch M "<<cylindervec[idx].cindex << " " <<
                      ccylvec[idx]->getCylinder()->_dcIndex << endl;
    }*/
//    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;
//    double timetaken = 0.0;
//    double time16 = 0.0;
//    double sqdisttermswithjustalpha;
//    vector<double> maxvec;
//    vector<double> minvec;
//    int rejects16 = 0;
//    int rejectsnavail =0;
//    mints = chrono::high_resolution_clock::now();
//    start
//    double bscoord[3* nbs * Ncyl];
//    double bsvec[nbs * Ncyl];
//        int totalneighbors = 0;
//    int counter1 = 0;
//    minte = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_vec(minte - mints);
//    std::cout<<"Vectorize time "<<elapsed_vec.count()<<endl;
    //V3 begins
//    int maxdistrejects = 0;
//    int IDrejects = 0;
//    _possibleBindingsstencil.clear();
//    mints = chrono::high_resolution_clock::now();
//    for(auto c : _compartment->getCylinders()) {
//
//        if(c->getType() != _filamentType) continue;
//
//        auto cc = c->getCCylinder();
//
//        for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
//            it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
//
//            //now re add valid binding sites
//            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
//                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
//
//                //loop through neighbors
//                //now re add valid based on CCNL
//                for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil
//                        (cc->getCylinder())) {
//
//                    if(cn->getParent() == c->getParent()) continue;
//                    if(cn->getType() != _filamentType) continue;
//
//                    auto ccn = cn->getCCylinder();
//
//                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
//                        it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
//
//                        if (areEqual(ccn->getCMonomer(*it2)->speciesBound(
//                                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0)) {
//
//                            //check distances..
//                            auto mp1 = (float)*it1 / SysParams::Geometry().cylinderNumMon[_filamentType];
//                            auto mp2 = (float)*it2 / SysParams::Geometry().cylinderNumMon[_filamentType];
//
//                            auto x1 = c->getFirstBead()->coordinate;
//                            auto x2 = c->getSecondBead()->coordinate;
//                            auto x3 = cn->getFirstBead()->coordinate;
//                            auto x4 = cn->getSecondBead()->coordinate;
//
//                            auto m1 = midPointCoordinate(x1, x2, mp1);
//                            auto m2 = midPointCoordinate(x3, x4, mp2);
//
//                            double dist = twoPointDistance(m1, m2);
//
//                            if(dist > _rMax || dist < _rMin) continue;
//
//                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
//                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);
//
//                            //add in correct order
//                            if(c->getID() > cn->getID())
//                                _possibleBindingsstencil.emplace(t1,t2);
//                        }
//                    }
//                }
//            }
//        }
//    }
//    minte = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_totalv3(minte - mints);
//    std::cout<<"Overall time vec v3 "<<elapsed_totalv3.count()<<endl;
//    std::cout<<"Tuple size "<<_possibleBindingsstencil.size()<<endl;
//    std::cout<<"maxdist rejects "<<maxdistrejects<<endl;
//    std::cout<<"ID rejects "<<IDrejects<<endl;
    //V3 ends

//    chrono::high_resolution_clock::time_point minIs, minIe, minroots, minroote ;
//    double ttup = 0.0;
//    double tsten = 0.0;
//    double tdense = 0.0;
    int counter2 =0; int counter3 = 0; int counter4 = 0; int counter5 = 0;
    int counter6 =0;int counter7 = 0; int counter8 = 0; int counter9 = 0;int counter10=0;
//    mints = chrono::high_resolution_clock::now();
    for(int i=0;i<Ncylincmp;i++){
        int cindex = cindexvec[i];
        cylinder c = cylindervec[cindex];
        if(c.type != _filamentType) {
            counter2++;
            continue;}
        double x1[3],x2[3];
        memcpy(x1, &coord[3*c.bindices[0]], 3 * sizeof(double));
        memcpy(x2, &coord[3*c.bindices[1]], 3 * sizeof(double));
        double X1X2[3] ={x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        //Check 2
/*
        for(auto cndummy:ncindices[i]){
            auto A =cylindervec[cndummy];
            if(ccylvec[A.cindex]->getCylinder()->getID() != A.ID)
            {
                std::cout<<"Mismatch in neighbors L of Cyl index "<<i<<" "
                        ""<<ccylvec[cndummy]->getCylinder()->getID()<<" "
                                 ""<<cylindervec[cndummy].ID<<endl;
            }
        }
*/

        int* cnindices = ncindices[i].data();
        for(int arraycount = 0; arraycount < ncindices[i].size();arraycount++){
            int cnindex = cnindices[arraycount];
//            std::cout<<*cnindex<<endl;
            cylinder cn = cylindervec[cnindex];
//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
// not contain ncs that will fail this cndn.
            if(c.filamentID == cn.filamentID){
                counter3++;
                continue;}
            if(c.type != cn.type){
                counter4++;
                 continue;}

            double x3[3], x4[3];
            memcpy(x3, &coord[3*cn.bindices[0]], 3 * sizeof(double));
            memcpy(x4, &coord[3*cn.bindices[1]], 3 * sizeof(double));
            double X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            double X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            double X1X3squared = sqmagnitude(X1X3);
            double X1X2squared = cylsqmagnitudevector[c.cindex];
            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            double X3X4squared = cylsqmagnitudevector[cn.cindex];
            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            for(int pos1 =0; pos1<nbs;pos1++) {
            //now re add valid binding sites
                if (areEqual(boundstate[2][offset + maxnbs *c.cindex + pos1], 1.0)) {
                auto mp1 = bindingsites.at(pos1);
                double A = X3X4squared;
                double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 * X1X3dotX1X2;
                double C1 = C - _rMinsq;
                double C2 = C - _rMaxsq;
                double b2m4ac1 = B*B - 4*A*C1;
                double b2m4ac2 = B*B - 4*A*C2;
                status1 = b2m4ac1 < 0;
                status2 = b2m4ac2 < 0;
                if(status1 && status2) {
                    counter6++;
                    continue;}
                if(!status1){
                    min1 = (-B + sqrt(b2m4ac1))/(2*A);
                    min2 = (-B - sqrt(b2m4ac1))/(2*A);
                    if(min1<min2) {
                        minveca[0] = (min1);
                        minveca[1] = (min2);
                    }
                    else{
                        minveca[0] = (min2);
                        minveca[1] = (min1);
                    }
                    if(minveca[0]< minparamcyl2 && minveca[1] > maxparamcyl2) {
                        counter7++;
                        continue;
                    }
                }
                if(!status2){
                    max1 = (-B + sqrt(b2m4ac2))/(2*A);
                    max2 = (-B - sqrt(b2m4ac2))/(2*A);
                    if(max1<max2) {
                        maxveca[0] = (max1);
                        maxveca[1] = (max2);
                    }
                    else{
                        maxveca[0] = (max2);
                        maxveca[1] = (max1);
                    }
                    if(maxveca[0] > maxparamcyl2 || maxveca[1] <minparamcyl2){
                        counter8++;
                        continue;
                    }
                }
                    for(int pos2 = 0; pos2<nbs;pos2++){
                    if (areEqual(boundstate[2][offset + maxnbs *cn.cindex + pos2], 1.0)) {
                        total++;
                        //check distances..
                        auto mp2 = bindingsites.at(pos2);
                        if(!status2) {
                            if (mp2 < maxveca[0] || mp2 > maxveca[1]) {
                                {
                                    counter9++;
                                continue;}
                            }
                        }
                        if(!status1){
                            if (mp2 > minveca[0] && mp2 < minveca[1]) {
                                {
                                    counter10++;
                                    continue;}
                            }
                        }
                        accepts++;

                        auto it1 = SysParams::Chemistry().bindingSites[_filamentType][pos1];
                        auto it2 = SysParams::Chemistry().bindingSites[_filamentType][pos2];
                        /*auto xx1 = ccylvec[cindex]->getCylinder();
                        auto xx2 = ccylvec[cnindex]->getCylinder();
//                        std::cout<<xx1->getID()<<" "<<c.ID<<" "<<xx2->getID()<<" "
//                                ""<<cn.ID<<endl;
                        if(xx1->getFirstBead()->_dbIndex != c.bindices[0] ||
                                xx1->getSecondBead()->_dbIndex != c.bindices[1])
                            std::cout<<"DB1 "<<xx1->getFirstBead()->_dbIndex<<" "
                                ""<<xx1->getSecondBead()->_dbIndex<<" "<<c
                                         .bindices[0]<<" "<<c.bindices[1]<<endl;
                        if(xx2->getFirstBead()->_dbIndex != cn.bindices[0] ||
                           xx2->getSecondBead()->_dbIndex != cn.bindices[1])
                        std::cout<<"DB2 "<<xx2->getFirstBead()->_dbIndex<<" "
                                ""<<xx2->getSecondBead()->_dbIndex<<" "<<cn
                                         .bindices[0]<<" "<<cn.bindices[1]<<endl;*/
                        auto t1 = tuple<CCylinder *, short>(ccylvec[cindex], it1);
                        auto t2 = tuple<CCylinder *, short>(ccylvec[cnindex], it2);
                        //add in correct order
                            _possibleBindingsstencil.emplace(t1, t2);
                    }
                }
            }
            }
        }
    }
//    minte = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_total4(minte - mints);
//    std::cout<<"Overall time vec "<<elapsed_total4.count()<<endl;
//    std::cout<<"Tuple time "<<ttup<<endl;
//    std::cout<<"Stencil time "<<tsten<<endl;
//    std::cout<<"Dense time "<<tdense<<endl;
//    std::cout<<"Tuple size "<<_possibleBindingsstencil.size()<<endl;
/*    std::cout<<"Total binding site pairs "<<16*totalneighbors<<endl;
    std::cout<<"Count cID < cnID "<<16*counter1<<endl;
    std::cout<<"Count ftype "<<16*counter2<<endl;
    std::cout<<"Count fID "<<16*counter3<<endl;
    std::cout<<"Count Ctype "<<16*counter4<<endl;
    std::cout<<"Count maxdist "<<16*counter5<<endl;
    std::cout<<"Count b2-4ac "<<4*counter6<<endl;
    std::cout<<"Count min out of range "<<4*counter7<<endl;
    std::cout<<"Count max out of range "<<4*counter8<<endl;
    std::cout<<"Count bs2 min out of range "<<counter9<<endl;
    std::cout<<"Count bs2 max out of range "<<counter10<<endl;
    std::cout<<"accepts "<<accepts<<endl;*/

//    accepts = 0;
//    _possibleBindingsstencil.clear();
//    mints = chrono::high_resolution_clock::now();
//    for(auto c : _compartment->getCylinders()) {
//        if (c->getType() != _filamentType) continue;
//
//        auto x1 = c->getFirstBead()->coordinate;
//        auto x2 = c->getSecondBead()->coordinate;
//        auto cc = c->getCCylinder();
//        vector<double> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
//
//        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(c)) {
//
//            if(cn->getParent() == c->getParent()) continue;
//            if(cn->getType() != _filamentType) continue;
//            if(c->getID() < cn->getID()) continue;
//            auto ccn = cn->getCCylinder();
//            auto x3 = cn->getFirstBead()->coordinate;
//            auto x4 = cn->getSecondBead()->coordinate;
//
//            vector<double> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
//            vector<double> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
//            double maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);
//
//            double mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
//                                                                                    x3,x4)));
//            mindistsq = mindistsq * mindistsq;
//            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;
//
//            double X1X3squared = sqmagnitude(X1X3);
//            double X1X2squared = cylsqmagnitudevector.at(c->_dcIndex);
//            double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
//            double X3X4squared = cylsqmagnitudevector.at(cn->_dcIndex);
//            double X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
//            double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
////            mins2 = chrono::high_resolution_clock::now();
//            int i = -1;
//            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
//                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
//                i++;
//                //now re add valid binding sites
//                if (areEqual(boundstate[2][offset + SysParams::Chemistry()
//                             .bindingSites[_filamentType].size() *c->_dcIndex + i], 1.0)) {
//                    auto mp1 = bindingsites.at(i);
//                    double A = X3X4squared;
//                    double B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
//                    double C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 * X1X3dotX1X2;
//                    double C1 = C - _rMinsq;
//                    double C2 = C - _rMaxsq;
//                    double b2m4ac1 = B*B - 4*A*C1;
//                    double b2m4ac2 = B*B - 4*A*C2;
//                    status1 = b2m4ac1 < 0;
//                    status2 = b2m4ac2 < 0;
//                    if(status1 && status2) continue;
//                    maxvec.clear();
//                    minvec.clear();
//                    if(!status1){
//                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
//                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
//                        if(min1<min2) {
//                            minvec.push_back(min1);
//                            minvec.push_back(min2);
//                        }
//                        else{
//                            minvec.push_back(min2);
//                            minvec.push_back(min1);
//                        }
//                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
//                            continue;
//                        }
//                    }
//                    if(!status2){
//                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
//                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
//                        if(max1<max2) {
//                            maxvec.push_back(max1);
//                            maxvec.push_back(max2);
//                        }
//                        else{
//                            maxvec.push_back(max2);
//                            maxvec.push_back(max1);
//                        }
//                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
//                            continue;
//                        }
//                    }
//                    int j =-1;
//                    for(auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
//                        it2 != SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {
//                        j++;
//                        bool check2 = true;
//                        if (areEqual(boundstate[2][offset + SysParams::Chemistry()
//                                                                    .bindingSites[_filamentType]
//                                                                    .size()*cn->_dcIndex + j], 1.0)) {
//                            total++;
//                            //check distances..
//                            auto mp2 = bindingsites.at(j);
//
//                            if(!status2) {
//                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1)) {
//                                    continue;
//                                }
//                            }
//                            if(!status1){
//                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1)) {
//                                    continue;
//                                }
//                            }
//
//                            accepts++;
//
//                            if(check2) {
//
//                                auto t1 = tuple<CCylinder *, short>(cc, *it1);
//                                auto t2 = tuple<CCylinder *, short>(ccn, *it2);
//
//                                //add in correct order
//                                if (c->getID() > cn->getID()) {
//                                    _possibleBindingsstencil.emplace(t1, t2);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//        }
//    }
//    //end
//    minte = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_total3(minte - mints);
//    std::cout<<"Overall time "<<elapsed_total3.count()<<endl;
//    std::cout<<"Total "<<total<<" accepts "<<accepts<<endl;
//    std::cout<<"Tuple size "<<_possibleBindingsstencil.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();
    updateBindingReaction(oldN, newN);
    /*std::cout<<"Motor consistency "<<isConsistent()<<endl;*/
}
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        removePossibleBindingsstencil(cc, *bit);
}
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#ifdef HYBRID_NLSTENCILLIST
//    cout<<"Removing "<<cc->getCylinder()->getID()<<" "<<bindingSite<<endl;
    auto HManager = _compartment->getHybridBindingSearchManager();
        HManager->removePossibleBindingsstencil(_idvec, cc, bindingSite);
/*    for(auto C:SubSystem::getstaticgrid()->getCompartments()){
        C->getHybridBindingSearchManager()->checkoccupancySIMD(_idvec);
    }*/
#else

    if(cc->getType() != _filamentType) return;

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencil.erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencil.begin(); it != _possibleBindingsstencil.end();) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
        _possibleBindingsstencil.erase(it++);

        else ++it;
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);


    //remove all neighbors which have this binding site pair
/*#ifdef HYBRID_NLSTENCILLIST
    auto C = cc->getCompartment();
    for (auto cn : C->getHybridBindingSearchManager()->getHNeighbors(cc->getCylinder(), HNLID)){
        if(cn->getType() != _filamentType) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (MotorBindingManager*)cn->getCompartment()->
            getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
            affectedManagers.push_back(m);
        }
    }
#else*/
    for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

        if(cn->getType() != _filamentType) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (MotorBindingManager*)cn->getCompartment()->
                    getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
//#endif
    //remove, update affected
    for(auto m : affectedManagers) {

        for (auto it = m->_possibleBindingsstencil.begin(); it !=
                m->_possibleBindingsstencil.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindingsstencil.erase(it++);

            else ++it;
        }

        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }
#endif
}
void MotorBindingManager::crosscheck(){
    cout<<"Motor NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
        <<_possibleBindingsstencil.size()<<endl;
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
    std::cout<<"Motor possible bindings size "<<_possibleBindings.size()<<" Total matches"
            " "<<
             matches<<endl;
    if(_possibleBindings.size() != matches || _possibleBindings.size() !=
       _possibleBindingsstencil.size()){
        cout<<"Motor. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
vector<tuple<CCylinder*, short>> MotorBindingManager::chooseBindingSitesstencil() {
#ifdef HYBRID_NLSTENCILLIST
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->chooseBindingSitesstencil(_idvec);
#else
    assert((_possibleBindingsstencil.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindingsstencil.size() - 1);
    auto it = _possibleBindingsstencil.begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}
#endif
#ifdef CUDAACCL_NL
void MotorBindingManager::assigncudavars() {

//    if(gpu_numpairs == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)), "cuda data transfer", "BindingManager.cu");
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
SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CylinderCylinderNL*> LinkerBindingManager::_neighborLists;
vector<CylinderCylinderNL*> MotorBindingManager::_neighborLists;
short LinkerBindingManager::HNLID;
short MotorBindingManager::HNLID;
short LinkerBindingManager::_idvec[2];
short MotorBindingManager::_idvec[2];
