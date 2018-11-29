
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
#ifdef HYBRID_NLSTENCILLIST
#include "Compartment.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "HybridBindingSearchManager.h"
#include "HybridNeighborListImpl.h"
#include "MotorGhost.h"
#include "MathFunctions.h"
#include "GController.h"
#include "SysParams.h"
#include "CUDAcommon.h"

vector<short> HybridBindingSearchManager::HNLIDvec;
using namespace mathfunc;
HybridBindingSearchManager::HybridBindingSearchManager(Compartment* compartment){
    _compartment = compartment;
    totaluniquefIDpairs = 0;
};
void HybridBindingSearchManager::setbindingsearchparameter
        (FilamentBindingManager* fmanager, short bstatepos, short ftype1, short
                        ftype2, float rMax, float rMin){
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> temp;
    bool isfound = false;
    vector<short> ftypepairs;
    if(ftype1 < ftype2)
        ftypepairs = {ftype1,ftype2};
    else
        ftypepairs = {ftype2, ftype1};
    for(short idx = 0; idx < totaluniquefIDpairs; idx++){
        vector<short> fIDpair = _filamentIDvec[idx];
        if(fIDpair[0] != ftypepairs[0] || fIDpair[1] != ftypepairs[1]) continue;
        else {
            isfound = true;
            _rMaxsqvec[idx].push_back(rMax *rMax);
            _rMinsqvec[idx].push_back(rMin *rMin);
            fManagervec[idx].push_back(fmanager);
            _possibleBindingsstencilvec[idx].push_back(temp);
            bstateposvec[idx].push_back(bstatepos);
            break;
        }
    }
    if(isfound == false){
        vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*,
                short>>> temp2;
        temp2.push_back(temp);
        vector<float> localrmaxsq ={rMax * rMax};
        vector<float> localrminsq = {rMin * rMin};
        vector<FilamentBindingManager*> localfmanager;
        vector<short> localbstateposvec = {bstatepos};
        _rMaxsqvec.push_back(localrmaxsq);
        _rMinsqvec.push_back(localrminsq);
        localfmanager.push_back(fmanager);
        fManagervec.push_back(localfmanager);
        _filamentIDvec.push_back(ftypepairs);
        _possibleBindingsstencilvec.push_back(temp2);
        bstateposvec.push_back(localbstateposvec);
        vector<double> bs1, bs2;
        vector<float> minvec = {(float)*(SysParams::Chemistry().bindingSites[ftypepairs[0]]
                                .begin())/ SysParams::Geometry().cylinderNumMon[ftypepairs[0]],
                                (float)*(SysParams::Chemistry().bindingSites[ftypepairs[1]]
                                .begin())/ SysParams::Geometry().cylinderNumMon[ftypepairs[1]]};
        vector<float> maxvec = {(float)*(SysParams::Chemistry().bindingSites[ftypepairs[0]]
                                .end() -1)/ SysParams::Geometry().cylinderNumMon[ftypepairs[0]],
                                (float)*(SysParams::Chemistry().bindingSites[ftypepairs[1]]
                                 .end() -1)/ SysParams::Geometry().cylinderNumMon[ftypepairs[1]]};
        for(auto it1 = SysParams::Chemistry().bindingSites[ftypepairs[0]].begin();
            it1 != SysParams::Chemistry().bindingSites[ftypepairs[0]].end(); it1++) {
            bs1.push_back((float)*it1 / SysParams::Geometry().cylinderNumMon[ftypepairs[0]]);
        }
        for(auto it1 = SysParams::Chemistry().bindingSites[ftypepairs[1]].begin();
            it1 != SysParams::Chemistry().bindingSites[ftypepairs[1]].end(); it1++) {
            bs2.push_back((float)*it1 / SysParams::Geometry().cylinderNumMon[ftypepairs[1]]);
        }
        minparamcyl2.push_back(minvec);
        maxparamcyl2.push_back(maxvec);
        bindingsites1.push_back(bs1);
        bindingsites2.push_back(bs1);
        totaluniquefIDpairs++;
    }
}
void HybridBindingSearchManager::addPossibleBindingsstencil(short idvec[2], CCylinder* cc,
                                    short bindingSite) {
//    std::cout<<"Adding "<<cc->getCylinder()->getID()<<" "<<bindingSite<<endl;
    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fIDpair = _filamentIDvec[idx].data();
    short bstatepos = bstateposvec[idx][idx2];
    short _filamentType = cc->getType();
    short HNLID = HNLIDvec[idx];
    short complimentaryfID;
    if(_filamentType != fIDpair[0] || _filamentType != fIDpair[1] ) return;
    else if(_filamentType == fIDpair[0]) complimentaryfID = fIDpair[1];
    else complimentaryfID = fIDpair[0];
    bool bstatecheck = false;
    float _rMaxsq = _rMaxsqvec[idx][idx2];
    float _rMinsq = _rMinsqvec[idx][idx2];
    if(bstatepos == 1)
        bstatecheck = areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);
    if(bstatepos == 2)
        bstatecheck = areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);

    //if we change other managers copy number
    vector<HybridBindingSearchManager*> affectedHManagers;
//    vector<MotorBindingManager*> affectedManagersM;
//    vector<LinkerBindingManager*> affectedManagersL;

    //add valid binding sites
    if (bstatecheck) {
        //loop through neighbors
        //now re add valid based on CCNL
        vector<Cylinder*> Neighbors;
        Neighbors = _HneighborList->getNeighborsstencil(HNLID, cc->getCylinder());
        for (auto cn : Neighbors) {
            Cylinder* c = cc->getCylinder();
            short _nfilamentType = cn->getType();
            if(_nfilamentType != complimentaryfID) return;
            if(cn->getParent() == c->getParent()) continue;

            auto ccn = cn->getCCylinder();

            for(auto it = SysParams::Chemistry().bindingSites[_nfilamentType].begin();
                it != SysParams::Chemistry().bindingSites[_nfilamentType].end(); it++) {
                bool bstatecheckn = false;
                if(bstatepos == 1)
                    bstatecheckn = areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_nfilamentType])->getN(), 1.0);
                if(bstatepos == 2)
                    bstatecheckn = areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_nfilamentType])->getN(), 1.0);

                if (bstatecheckn) {

                    //check distances..
                    auto mp1 = (float)bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
                    auto mp2 = (float)*it / SysParams::Geometry().cylinderNumMon[_nfilamentType];

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
                        _possibleBindingsstencilvec[idx][idx2].emplace(t1,t2);
                    }

                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {

                            _possibleBindingsstencilvec[idx][idx2].emplace(t2,t1);
                        }
                            //add in other
                        else {

                            auto m = cn->getCompartment()->getHybridBindingSearchManager();

                            affectedHManagers.push_back(m);

                            m->_possibleBindingsstencilvec[idx][idx2].emplace(t2,t1);
                        }
                    }
                }
            }
        }
    }

    //update affected
    for(auto m : affectedHManagers) {
        int newNOther = m->_possibleBindingsstencilvec[idx][idx2].size();
        m->fManagervec[idx][idx2]->updateBindingReaction(newNOther);
    }

    //update this manager
    int newN = _possibleBindingsstencilvec[idx][idx2].size();
    fManagervec[idx][idx2]->updateBindingReaction(newN);
}
void HybridBindingSearchManager::removePossibleBindingsstencil(short idvec[2], CCylinder*
                                    cc, short bindingSite) {
    /*std::cout<<"Removing possible bindings of cylinder with ID "<<cc->getCylinder()
            ->getID()<<" with cindex "<<cc->getCylinder()->_dcIndex<<" "<<endl;*/
    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fIDpair = _filamentIDvec[idx].data();
    short bstatepos = bstateposvec[idx][idx2];
    short _filamentType = cc->getType();
    short HNLID = HNLIDvec[idx];

    if(_filamentType != fIDpair[0] && _filamentType != fIDpair[1] ) return;

    //if we change other managers copy number
    vector<HybridBindingSearchManager*> affectedHManagers;
//    int oldN = _possibleBindingsstencilvec[idx][idx2].size();
    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencilvec[idx][idx2].erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencilvec[idx][idx2].begin(); it !=
            _possibleBindingsstencilvec[idx][idx2].end();) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {

            _possibleBindingsstencilvec[idx][idx2].erase(it++);
        }

        else ++it;
    }

    int newN = _possibleBindingsstencilvec[idx][idx2].size();
    /*std::cout<<"Cyl "<<cc->getCylinder()->getID()<<" bs "<<bindingSite<<endl;
    std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<_compartment->coordinates()
    [1]<<" "<<_compartment->coordinates()[2]<<" pairs removed "<<oldN-newN<<endl;*/

    fManagervec[idx][idx2]->updateBindingReaction(newN);
    //remove all neighbors which have this binding site pair
/*    auto Neighbors = _HneighborList->getNeighborsstencil(HNLID, cc->getCylinder());
    for (auto cn : Neighbors){
        short _nfilamentType = cn->getType();
        if(_nfilamentType != fIDpair[0] && _nfilamentType != fIDpair[1] ) return;
        //inaccurate. replace later.
        if(cn->getCompartment() != _compartment) {

            auto m = cn->getCompartment()->getHybridBindingSearchManager();

            if(find(affectedHManagers.begin(), affectedHManagers.end(), m) ==
                    affectedHManagers.end())
                affectedHManagers.push_back(m);
        }
    }

    //remove, update affected
    for(auto m : affectedHManagers) {*/
    //Go through enclosing compartments to remove all entries with the current cylinder
    // and binding site as values.
    for(auto nc: _compartment->getenclosingNeighbours()){
        if(nc != _compartment) {
            auto m = nc->getHybridBindingSearchManager();
//        int oldN = m->_possibleBindingsstencilvec[idx][idx2].size();
            for (auto it = m->_possibleBindingsstencilvec[idx][idx2].begin(); it !=
                 m->_possibleBindingsstencilvec[idx][idx2].end();) {

                if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                    m->_possibleBindingsstencilvec[idx][idx2].erase(it++);

                else ++it;
            }
            int newNOther = m->_possibleBindingsstencilvec[idx][idx2].size();
            m->fManagervec[idx][idx2]->updateBindingReaction(newNOther);
            /*std::cout<<"Cmp "<<m->_compartment->coordinates()[0]<<" "
                    ""<<m->_compartment->coordinates()
            [1]<<" "<<m->_compartment->coordinates()[2]<<" pairs removed "
                    ""<<oldN-newNOther<<endl;*/
        }
    }
}

void HybridBindingSearchManager::checkoccupancy(short idvec[2]){
    short idx = idvec[0];
    short idx2 = idvec[1];
    auto pbs = _possibleBindingsstencilvec[idx][idx2];
    for(auto pair = pbs.begin(); pair != pbs.end(); pair++){
        auto leg1 = pair->first;
        auto leg2 = pair->second;
        CCylinder* ccyl1 = get<0>(leg1);
        short bs1 = get<1>(leg1);
        CCylinder* ccyl2 = get<0>(leg2);
        short bs2 = get<1>(leg2);
        bool state1 = areEqual(ccyl1->getCMonomer(bs1)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[0])->getN(), 1.0);
        bool state2 = areEqual(ccyl2->getCMonomer(bs2)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[0])->getN(), 1.0);
        if(state1 != true || state2 != true){
            std::cout<<"OOPS occupied species exist "<<state1<<" "<<state2<<endl;
            SpeciesBound* sm1 = ccyl1->getCMonomer(bs1)->speciesMotor(0);
            SpeciesBound* sm2 = ccyl2->getCMonomer(bs2)->speciesMotor(0);
            SpeciesBound* BM1 = ccyl1->getCMonomer(bs1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[0]);
            SpeciesBound* BM2 = ccyl2->getCMonomer(bs2)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[0]);
            std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<_compartment->coordinates()
            [1]<<" "<<_compartment->coordinates()[2]<<" Motor "<<ccyl1->getCylinder()->getID()<<" "<<bs1<<" "
                    ""<<ccyl2->getCylinder()->getID()<<" "<<
                     ""<<bs2<< endl;
            std::cout<<"Motor "<<sm1->getN()<<" "<<sm2->getN()<<" BOUND "<<BM1->getN()<<" "<<BM2->getN()<<endl;
//            exit(EXIT_FAILURE);
        }
    }
}
void HybridBindingSearchManager::updateAllPossibleBindingsstencil() {
    int counter1, counter2, counter3, counter4, counter5, counter6, counter7, counter8,
            counter9, counter10;
    counter1 = 0;counter2 = 0;counter3 = 0;counter4 = 0;counter5 = 0;counter6 = 0;
    counter7 = 0;counter8 = 0;counter9 = 0;counter10 = 0;
    int totalneighbors = 0;
    for (int idx = 0; idx < totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
    for (int idx2 = 0; idx2 < countbounds; idx2++) {
        _possibleBindingsstencilvec[idx][idx2].clear();
    }
    }

    double min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    double minveca[2];
    double maxveca[2];
    int accepts = 0;
    int total = 0;
    int Ncyl = Cylinder::getCylinders().size();
    double* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    double *coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int Ncylincmp = _compartment->getCylinders().size();
    int cindexvec[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int> ncindex; //helper vector
    int maxnbs = SysParams::Mechanics().maxbindingsitespercylinder;

    int idx; int idx2;
    for(idx = 0; idx<totaluniquefIDpairs;idx++) {
        long id = 0;
        ncindices.clear();
        auto fpairs = _filamentIDvec[idx].data();
        int offset1 = SysParams::Mechanics().bsoffsetvec.at(fpairs[0]);
        int offset2 = SysParams::Mechanics().bsoffsetvec.at(fpairs[1]);
        int nbs1 = SysParams::Chemistry().bindingSites[fpairs[0]].size();
        int nbs2 = SysParams::Chemistry().bindingSites[fpairs[1]].size();

        int Nbounds = _rMaxsqvec[idx].size();
//        auto bindingsitevec1 = SysParams::Chemistry().bindingSites[fpairs[0]];
//        auto bindingsitevec2 = SysParams::Chemistry().bindingSites[fpairs[1]];
        for (auto c : _compartment->getCylinders()) {
            cindexvec[id] = c->_dcIndex;
//        totalneighbors += _neighborLists[_nlIndex]->getNeighborsstencil(c).size();
            id++;
            auto Neighbors = _HneighborList->getNeighborsstencil(HNLIDvec[idx], c);
            ncindex.reserve(Neighbors.size());
            for (auto cn : Neighbors) {
                if (c->getID() > cn->getID())
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
                      ccylvec[a1]->getCylinder()->_dcIndex << endl;
    }*/
        for (int i = 0; i < Ncylincmp; i++) {
            int cindex = cindexvec[i];
            cylinder c = cylindervec[cindex];
            short complimentaryfID;
            if (c.type != fpairs[0] && c.type != fpairs[1]) {
            counter2++;
                continue;
            }
            else if(c.type == fpairs[0]) complimentaryfID = fpairs[1];
            else complimentaryfID = fpairs[0];
            double x1[3], x2[3];
            memcpy(x1, &coord[3 * c.bindices[0]], 3 * sizeof(double));
            memcpy(x2, &coord[3 * c.bindices[1]], 3 * sizeof(double));
            double X1X2[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

            int *cnindices = ncindices[i].data();
            for (int arraycount = 0; arraycount < ncindices[i].size(); arraycount++) {
                int cnindex = cnindices[arraycount];
                cylinder cn = cylindervec[cnindex];
//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
//              not contain ncs that will fail this cndn.
                if (c.filamentID == cn.filamentID) {
                counter3++;
                    continue;
                }
                if(c.type != complimentaryfID) continue;
//                if (c.type != cn.type) {
////                counter4++;
//                    continue;
//                }

                double x3[3], x4[3];
                memcpy(x3, &coord[3 * cn.bindices[0]], 3 * sizeof(double));
                memcpy(x4, &coord[3 * cn.bindices[1]], 3 * sizeof(double));
                double X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
                double X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
                double X1X3squared = sqmagnitude(X1X3);
                double X1X2squared = cylsqmagnitudevector[c.cindex];
                double X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
                double X3X4squared = cylsqmagnitudevector[cn.cindex];
                double X1X3dotX3X4 = scalarprojection(X1X3, X3X4);
                double X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
                //Number of binding sites on the cylinder/filament of Type A.
                for (int pos1 = 0; pos1 < nbs1; pos1++) {
                    //Number of binding distance pairs we are looking at
                    for (idx2 = 0; idx2 < Nbounds; idx2++) {
                        short bstatepos = bstateposvec[idx][idx2];
                    //now re add valid binding sites
                    if (areEqual(boundstate[bstatepos][offset1 + maxnbs * c.cindex + pos1], 1.0)) {
                        auto mp1 = bindingsites1[idx][pos1];
                        double A = X3X4squared;
                        double B = 2.0 * X1X3dotX3X4 - 2.0 * mp1 * X3X4dotX1X2;
                        double C = X1X3squared + mp1 * mp1 * X1X2squared -
                                   2.0 * mp1 * X1X3dotX1X2;
                        double Bsq = B*B;

                            double C1 = C - _rMinsqvec[idx][idx2];
                            double C2 = C - _rMaxsqvec[idx][idx2];
                            double b2m4ac1 = Bsq - 4 * A * C1;
                            double b2m4ac2 = Bsq - 4 * A * C2;
                            status1 = b2m4ac1 < 0;
                            status2 = b2m4ac2 < 0;
                            if (status1 && status2) {
                              counter6++;
                                continue;
                            }
                            if (!status1) {
                                min1 = (-B + sqrt(b2m4ac1)) / (2 * A);
                                min2 = (-B - sqrt(b2m4ac1)) / (2 * A);
                                if (min1 < min2) {
                                    minveca[0] = (min1);
                                    minveca[1] = (min2);
                                } else {
                                    minveca[0] = (min2);
                                    minveca[1] = (min1);
                                }
                                //Compare the MIN solutions are within the first and last
                                // binding sites in filament/cylinder of type B.
                                if (minveca[0] < minparamcyl2[idx][1] &&
                                    minveca[1] > maxparamcyl2[idx][1]) {
                                  counter7++;
                                    continue;
                                }
                            }
                            if (!status2) {
                                max1 = (-B + sqrt(b2m4ac2)) / (2 * A);
                                max2 = (-B - sqrt(b2m4ac2)) / (2 * A);
                                if (max1 < max2) {
                                    maxveca[0] = (max1);
                                    maxveca[1] = (max2);
                                } else {
                                    maxveca[0] = (max2);
                                    maxveca[1] = (max1);
                                }
                                //Compare the mAX solutions are within the first and last
                                // binding sites in filament/cylinder of type B.
                                if (maxveca[0] > maxparamcyl2[idx][1] ||
                                    maxveca[1] < minparamcyl2[idx][1]) {
                                  counter8++;
                                    continue;
                                }
                            }
                            for (int pos2 = 0; pos2 < nbs2; pos2++) {
                                if (areEqual(boundstate[bstatepos][offset2 + maxnbs * cn.cindex + pos2], 1.0)) {
                                    total++;
                                    //check distances..
                                    auto mp2 = bindingsites2[idx][pos2];
                                    if (!status2) {
                                        if (mp2 < maxveca[0] || mp2 > maxveca[1]) {
                                            {
                                    counter9++;
                                                continue;
                                            }
                                        }
                                    }
                                    if (!status1) {
                                        if (mp2 > minveca[0] && mp2 < minveca[1]) {
                                            {
                                    counter10++;
                                                continue;
                                            }
                                        }
                                    }
                                    accepts++;
                                    auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][pos1];
                                    auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][pos2];
                                    //Checker
                                    /*bool state1 = false;bool state2 = false;
                                    if(bstatepos == 2) {
                                        state1 = areEqual(ccylvec[cindex]->getCMonomer(
                                                it1)->speciesBound(
                                                SysParams::Chemistry().motorBoundIndex[0])->getN(),
                                                          1.0);
                                        state2 = areEqual(ccylvec[cnindex]->getCMonomer(
                                                it2)->speciesBound(
                                                SysParams::Chemistry().motorBoundIndex[0])->getN(),
                                                          1.0);
                                    }
                                    else if(bstatepos == 1){
                                        state1 = areEqual(ccylvec[cindex]->getCMonomer(
                                                it1)->speciesBound(
                                                SysParams::Chemistry().linkerBoundIndex[0])
                                                                  ->getN(),
                                                          1.0);
                                        state2 = areEqual(ccylvec[cnindex]->getCMonomer(
                                                it2)->speciesBound(
                                                SysParams::Chemistry().linkerBoundIndex[0])
                                                                  ->getN(),
                                                          1.0);
                                    }
                                    if(state1 != true || state2 != true){
                                        std::cout<<"states "<<state1<<" "<<state2<<endl;
                                        std::cout<<"bstate " <<boundstate[bstatepos][offset1 + nbs1 *
                                                    c.cindex + pos1]<<" "
                                                ""<<boundstate[bstatepos][offset2 + nbs2
                                                    * cn.cindex + pos2]<<endl;
                                        std::cout<<"bstatepos "<<bstatepos<<endl;
                                        std::cout<<"Adding an occupied site "<<endl;
                                        std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<_compartment->coordinates()
                                        [1]<<" "<<_compartment->coordinates()[2]<<" Motor"
                                                " "<<ccylvec[cindex]->getCylinder()
                                                ->getID()<<" "<<it1<<" "
                                                         ""<<ccylvec[cnindex]->getCylinder()->getID()<<" "<<
                                                 ""<<it2<< endl;
                                    }*/
                                    auto t1 = tuple<CCylinder *, short>(ccylvec[cindex],
                                                                        it1);
                                    auto t2 = tuple<CCylinder *, short>(ccylvec[cnindex],
                                                                        it2);
                                    //add in correct order
                                    _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

//Place in the appropriate BindingManager
    for(idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (idx2 = 0; idx2 < countbounds; idx2++) {
            int newNOther = _possibleBindingsstencilvec[idx][idx2].size();
            fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }
    /*std::cout<<"Total binding site pairs "<<16*totalneighbors<<endl;
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

//    int oldN = _bindingSpecies->getN();
//    int newN = numBindingSitesstencil();
//    updateBindingReaction(oldN, newN);
}
void HybridBindingSearchManager::addtoHNeighborList(){
    HNLIDvec.reserve(totaluniquefIDpairs);
    for(short idx = 0; idx < totaluniquefIDpairs; idx++) {
        vector<float> temprMaxsq = _rMaxsqvec[idx];
        vector<float> temprMinsq = _rMinsqvec[idx];
        vector<short> ftypepair = _filamentIDvec[idx];
        float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftypepair[0]],
                                    SysParams::Geometry().cylinderSize[ftypepair[1]]);
        float maxrMaxsq = *max_element(temprMaxsq.begin(), temprMaxsq.end());
        float minrMinsq = *min_element(temprMinsq.begin(), temprMinsq.end());
        for(int idx2 = 0; idx2<temprMaxsq.size();idx2++) {
            short temp = _HneighborList->setneighborsearchparameters(ftypepair[0],
                                                                     ftypepair[1], false, true, localmaxcylsize +
                                                                                                sqrt(_rMaxsqvec[idx][idx2]),
                                                                     max( sqrt(_rMinsqvec[idx][idx2]) - localmaxcylsize,float(0.0)));
            HNLIDvec[idx]=(temp);
            short idvec[2] = {idx, idx2};
            fManagervec[idx][idx2]->setHNLID(HNLIDvec[idx], idvec);
//            std::cout << "HNLID " << HNLIDvec[idx]<< endl;
        }
    }
}

HybridCylinderCylinderNL* HybridBindingSearchManager::_HneighborList;
#endif