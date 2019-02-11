
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
/*#include <boost/range/counting_range.hpp>
#include <thrust/execution_policy.h>
#include <thrust/system/omp/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/iterator/iterator_traits.h>
#include <thrust/binary_search.h>
#include <cstdlib>
#include <thrust/scan.h>*/

vector<short> HybridBindingSearchManager::HNLIDvec;
using namespace mathfunc;
HybridBindingSearchManager::HybridBindingSearchManager(Compartment* compartment){
    _compartment = compartment;
    totaluniquefIDpairs = 0;
    googlepossiblel.set_empty_key(4294967295);
    googlepossiblel.set_deleted_key(4294967294);
    googlereversepossiblel.set_empty_key(4294967295);
    googlereversepossiblel.set_deleted_key(4294967294);
    //
    googlepossiblem.set_empty_key(4294967295);
    googlepossiblem.set_deleted_key(4294967294);
    googlereversepossiblem.set_empty_key(4294967295);
    googlereversepossiblem.set_deleted_key(4294967294);
    //
    googlepossibleORIG.set_empty_key(4294967295);
    googlepossibleORIG.set_deleted_key(4294967294);
    googlereversepossibleORIG.set_empty_key(4294967295);
    googlereversepossibleORIG.set_deleted_key(4294967294);
    uint N = 10000*5000;
    linker1.reserve(N);
    linker2.reserve(N);
    motor1.reserve(N);
    motor2.reserve(N);
    gdmap tempmap;
    for(int i =0 ; i < nthreads; i++) {
        vecmapL[2 * i] = (tempmap);
        vecmapL[2 * i + 1] = (tempmap);
        vecmapM[2 * i] = (tempmap);
        vecmapM[2 * i + 1] = (tempmap);
    }

    gdmap tempmapgd;
/*    for(int i =1 ; i <= 2 * nthreads; i++) {
        _gdmapstruct.getelement(i) = tempmapgd;
    }*/

    vector<uint32_t> temp;
    pairslinker = temp;
    pairsmotor = temp;
    /*pairslinker.push_back(temp);
    pairslinker.push_back(temp);
    pairsmotor.push_back(temp);
    pairsmotor.push_back(temp);*/

    vector<bool> temp2;
    pairvaluespecieslinker.push_back(temp2);
    pairvaluespecieslinker.push_back(temp2);
    pairvaluespeciesmotor.push_back(temp2);
    pairvaluespeciesmotor.push_back(temp2);



};

void HybridBindingSearchManager::setbindingsearchparameter
        (FilamentBindingManager* fmanager, short bstatepos, short ftype1, short
        ftype2, float rMax, float rMin){
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> temp;
    unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*, short>>> rtemp;
    unordered_map<uint32_t, vector<uint32_t>> tempuint;
    unordered_map<uint32_t, vector<uint32_t>> rtempuint;
    gdmap tempgoogle;
    gdmap rtempgoogle;
    int tempNbind;
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
            _reversepossibleBindingsstencilvec[idx].push_back(rtemp);
            _possibleBindingsstencilvecuint[idx].push_back(tempuint);
            _reversepossibleBindingsstencilvecuint[idx].push_back(rtempuint);
            googlepossible[idx].push_back(tempgoogle);
            googlereversepossible[idx].push_back(rtempgoogle);
            bstateposvec[idx].push_back(bstatepos);
            Nbindingpairs[idx].push_back(tempNbind);
            //Set deleted key and initial key
            googlepossible[idx][googlepossible[idx].size()-1].set_empty_key(4294967295);
            googlepossible[idx][googlepossible[idx].size()-1].set_deleted_key(4294967294);
            googlereversepossible[idx][googlereversepossible[idx].size()-1].set_empty_key
            (4294967295);
            googlereversepossible[idx][googlereversepossible[idx].size()-1].set_deleted_key
            (4294967294);

            break;
        }
    }
    if(isfound == false){
        vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*,
                short>>> temp2;
        vector<unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*, short>>>>
        rtemp2;
        vector<unordered_map<uint32_t, vector<uint32_t>>> temp2uint;
        vector<unordered_map<uint32_t, vector<uint32_t>>> rtemp2uint;
        vector<gdmap> temp2google;
        vector<gdmap> rtemp2google;
        vector<int> tempNbind2;
        temp2.push_back(temp);
        rtemp2.push_back(rtemp);
        temp2uint.push_back(tempuint);
        rtemp2uint.push_back(rtempuint);
        temp2google.push_back(tempgoogle);
        rtemp2google.push_back(rtempgoogle);
        Nbindingpairs.push_back(tempNbind2);
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
        _reversepossibleBindingsstencilvec.push_back(rtemp2);
        _possibleBindingsstencilvecuint.push_back(temp2uint);
        _reversepossibleBindingsstencilvecuint.push_back(rtemp2uint);
        googlepossible.push_back(temp2google);
        googlereversepossible.push_back(rtemp2google);

        // set deleted key and insert key
        googlepossible[googlepossible.size() - 1][0].set_empty_key(4294967295);
        googlepossible[googlepossible.size() - 1][0].set_deleted_key(4294967294);
        googlereversepossible[googlereversepossible.size()-1][0].set_empty_key(4294967295);
        googlereversepossible[googlereversepossible.size()-1][0].set_deleted_key(4294967294);

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
#ifdef SIMDBINDINGSEARCH
    if(SysParams::INITIALIZEDSTATUS) {
        bool LinkerorMotor = true;
        if (idvec[1] == 1)
            LinkerorMotor = false;

        for (auto ncmp:_compartment->getenclosingNeighbours()) {
            if (LinkerorMotor)
                ncmp->SIMDcoordinates4linkersearch(0);
            else
                ncmp->SIMDcoordinates4motorsearch(0);
        }
        //Generate current coordinate
        int N = 8;
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(
                N), bindsitecoordinatesZ(N);
        vector<uint16_t> cindex_bs(N);
        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        auto cyl = cc->getCylinder();
        uint cindex = cyl->_dcIndex;
        auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];

        if (checkftype)
            _filamentType = cylinderstruct.type;
        auto x1 = cyl->getFirstBead()->coordinate;
        auto x2 = cyl->getSecondBead()->coordinate;
        auto mp = (float) bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];
        auto coord = midPointCoordinate(x1, x2, mp);
        for (uint k = 0; k < N; k++) {
            bindsitecoordinatesX[k] = coord[0];
            bindsitecoordinatesY[k] = coord[1];
            bindsitecoordinatesZ[k] = coord[2];
            //last 4 bits are binding site while first 12 bits are cylinder index.
            cindex_bs[k] = 1;
        }
        //Create input vector for SIMD calculations
        dist::Coords bscoords;
        bscoords.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                             bindsitecoordinatesZ, cindex_bs);

        //Calculate distances
        if (idvec[0] == 0 && idvec[1] == 0) {
            //Linker
            calculatebspairsforacylinder<true>(idvec, cindex, cc,
                                               bindingSite, bspairslinkerself,
                                               bspairslinker, bscoords);
        } else {
            //Motor
            calculatebspairsforacylinder<false>(idvec, cindex, cc,
                                                bindingSite, bspairsmotorself, bspairsmotor,
                                                bscoords);
        }

        //update affected
        short idx = idvec[0];
        short idx2 = idvec[1];
        //Calculate N
        minsfind = chrono::high_resolution_clock::now();
        for (short idx = 0; idx < totaluniquefIDpairs; idx++) {
            int countbounds = _rMaxsqvec[idx].size();
            for (short idx2 = 0; idx2 < countbounds; idx2++) {
                short idvec[2] = {idx, idx2};;
                countNpairsfound(idvec);
                fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);
            }
        }
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_countsites(minefind - minsfind);
        SIMDcountbs += elapsed_countsites.count();
    }
#else
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
                        _reversepossibleBindingsstencilvec[idx][idx2][t2].push_back(t1);
                    }
                    else {
                        //add in this compartment
                        if(cn->getCompartment() == _compartment) {
                            _possibleBindingsstencilvec[idx][idx2].emplace(t2,t1);
                            _reversepossibleBindingsstencilvec[idx][idx2][t1].push_back(t2);
                        }
                            //add in other
                        else {
                            auto m = cn->getCompartment()->getHybridBindingSearchManager();

                            affectedHManagers.push_back(m);

                            m->_possibleBindingsstencilvec[idx][idx2].emplace(t2,t1);
                            m->_reversepossibleBindingsstencilvec[idx][idx2][t1].push_back(t2);
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
#endif

}

template<bool LinkerorMotor>
void HybridBindingSearchManager::calculatebspairsforacylinder
        (short idvec[2], uint cIndex1, CCylinder* cCylinder, short bindingSite,
         dist::dOut<1, true>& bspairsoutSself,  dist::dOut<1, false>&
        bspairsoutS, dist::Coords& Coord) {
    const uint D = 1;
    short idx = idvec[0];
    short idx2 = idvec[1];


    //Get cylinder struct
    auto cyl = cCylinder->getCylinder();
    int cindex = cyl->_dcIndex;
    auto cylinderstruct = CUDAcommon::serlvars.cylindervec;
    auto cylinder1 = cylinderstruct[cindex];

    auto cylcmp1 = _compartment->Cyldcindexvec;
    short bstatepos = bstateposvec[idx][idx2];
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;

//    if(cc !=cCylinder){
//        cout<<"Huge mistake CCylinders don't match"<<endl;
//    }
    short _filamentType = cCylinder->getType();
    auto fIDpair = _filamentIDvec[idx].data();

    short complimentaryfID;
    if(_filamentType != fIDpair[0] || _filamentType != fIDpair[1] ) return;
    else if(_filamentType == fIDpair[0]) complimentaryfID = fIDpair[1];
    else complimentaryfID = fIDpair[0];

//    short position = 0;
//    for(short bs:SysParams::Chemistry().bindingSites[_filamentType]){
//        if(bs == bindingSite)
//            break;
//        else
//            position++;
//    }

    bool bstatecheck = false;
    if(bstatepos == 1)
        bstatecheck = areEqual(cCylinder->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);
    if(bstatepos == 2)
        bstatecheck = areEqual(cCylinder->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);

    //if the cylinder and bindingsite to be added are occupied, return.
    if(!bstatecheck) return;

    //This needs enclosing neighbors to ensure all neighbors are added
    for(auto ncmp: _compartment->getenclosingNeighbours()){
        //get bspairs based on idvec.
/*        cout<<"ncmpID "<<ncmp->getID()<<" simdSize "<<ncmp->getSIMDcoords<LinkerorMotor>()
                .size()<<" ncyls in ncmp "<<ncmp->getCylinders().size()<<" "<<LinkerorMotor<<endl;*/
            if(ncmp->getSIMDcoords<LinkerorMotor>().size()) {
                bspairsoutS.reset_counters();
                dist::find_distances(bspairsoutS, Coord,
                                     ncmp->getSIMDcoords<LinkerorMotor>(),
                                     t_serial);
                //Parse through dOut
                parsebspairsforacylinder<D, false>(idvec, bspairsoutS, ncmp, cindex,
                                                   bindingSite, complimentaryfID);
            }
    }
}

template <uint D, bool SELF>
void HybridBindingSearchManager::parsebspairsforacylinder(short idvec[2],
                                dist::dOut<D, SELF>&bspairsoutS, Compartment* ncmp,
                                uint cindex1, short bindingSite, short complimentaryfID){
    //Loop through output
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    auto cylinderstruct = CUDAcommon::serlvars.cylindervec;
    auto cylinder1 = cylinderstruct[cindex1];
    short idx = idvec[0];
    short idx2 = idvec[1];
    short bstatepos = bstateposvec[idx][idx2];
    //Trace data corresponding to cylinder1
    short _filamentType = cylinder1.type;
    auto it1 = bindingSite;
    auto pos = find (SysParams::Chemistry().bindingSites[_filamentType].begin(),
            SysParams::Chemistry().bindingSites[_filamentType].end(), it1);
    short j = pos - SysParams::Chemistry().bindingSites[_filamentType].begin();

    uint32_t shiftedCID = cylinder1.ID<<4;

    uint32_t t1 = shiftedCID|j;

    auto cylcmp2 = ncmp->Cyldcindexvec;
//    auto cylcmp1 = _compartment->Cyldcindexvec;
    short dim = D -1;
    uint N = bspairsoutS.counter[dim];

    for(uint pid = 0; pid < N; pid++) {

        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];
        uint16_t cmpcylidx2 = site2;
        cmpcylidx2 = (cmpcylidx2 >> 4);
        int cIndex2 = cylcmp2[cmpcylidx2];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];
        uint32_t cID2 = cylinder2.ID;
        uint32_t shiftedcID2 = cID2<<4;
        short bsite2 = mask & site2;

        //Check if other conditions are satisfied
        if (cylinder1.filamentID != cylinder2.filamentID &&
            (cylinder1.type == _filamentType && cylinder2.type == complimentaryfID)){
                //Add code to eliminate cylinders that are adjacent on the same filament.
            auto ccn = ccylvec[cIndex2];
            short _nfilamentType = ccn->getType();
            if(_nfilamentType != complimentaryfID) continue;
            auto it2 = SysParams::Chemistry().bindingSites[_nfilamentType][bsite2];

            bool bstatecheckn = false;
/*            if(ccn->getCylinder()->_dcIndex != cIndex2)
                cout<<"Neighbor Ccylinder does not represent the cylinder "
                      ""<<cIndex2<<" "<<ccn->getCylinder()->_dcIndex<<endl;*/
            if(bstatepos == 1)
                bstatecheckn = areEqual(ccn->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_nfilamentType])->getN(), 1.0);
            if(bstatepos == 2)
                bstatecheckn = areEqual(ccn->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_nfilamentType])->getN(), 1.0);


            if(bstatecheckn){
                uint32_t t2 = shiftedcID2|bsite2;
                googlepossible[idx][idx2][t1].push_back(t2);
                googlereversepossible[idx][idx2][t2].push_back(t1);
            }
            //}
        }
    }

}

void HybridBindingSearchManager::removePossibleBindingsstencil(short idvec[2], CCylinder*
                                    cc, short bindingSite) {
#ifdef SIMD_BINDINGSEARCH
    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fIDpair = _filamentIDvec[idx].data();
    short _filamentType = cc->getType();

    if(_filamentType != fIDpair[0] && _filamentType != fIDpair[1] ) return;

    //if we change other managers copy number
    vector<HybridBindingSearchManager*> affectedHManagers;

    //remove all tuples which have this ccylinder as key
    uint32_t t = cc->getCylinder()->getID()<<4;
    short pos = find (SysParams::Chemistry().bindingSites[_filamentType].begin(),
                      SysParams::Chemistry().bindingSites[_filamentType].end(),
                      bindingSite) -
                      SysParams::Chemistry().bindingSites[_filamentType].begin();
    t = t|pos;
    googlepossible[idx][idx2].erase(t);

    //remove all tuples which have this as value
    //Iterate through the reverse map
    auto keys = googlereversepossible[idx][idx2][t];//keys that contain t as
    // value in possiblebindings
    for(auto k:keys){
        //get the iterator range that corresponds to this key.
        auto range = googlepossible[idx][idx2].equal_range(k);
        //iterate through the range
        for(auto it = range.first; it != range.second;){
            //Go through the value vector and delete entries that match
            it->second.erase(remove(it->second.begin(), it->second.end(), t), it->second
            .end());
            it++;
        }
    }
    //remove from the reverse map.
    googlereversepossible[idx][idx2][t].clear();

    countNpairsfound(idvec);
    fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);

    //remove all neighbors which have this binding site pair
    //Go through enclosing compartments to remove all entries with the current cylinder
    // and binding site as values.
    for(auto nc: _compartment->getenclosingNeighbours()){
        if(nc != _compartment) {
            auto m = nc->getHybridBindingSearchManager();

            //Iterate through the reverse map
            auto keys = m->googlereversepossible[idx][idx2][t];//keys that
            // contain t as value in possiblebindings
            for(auto k:keys){
                //get the iterator range that corresponds to this key.
                auto range = m->googlepossible[idx][idx2].equal_range(k);
                //iterate through the range
                for(auto it = range.first; it != range.second;){
                    //Go through the value vector and delete entries that match
                    auto deleteiterator = remove(it->second.begin(),it->second.end(),t);
                    if(deleteiterator != it->second.end()){
                       it->second.erase(deleteiterator);
                    }
                    /*it->second.erase(remove(it->second.begin(), it->second.end(), t),
                            it->second.end());*/
                    it++;
                }
            }
            //remove from the reverse map.
            m->googlereversepossible[idx][idx2][t].clear();

            m->countNpairsfound(idvec);
            m->fManagervec[idx][idx2]->updateBindingReaction(m->Nbindingpairs[idx][idx2]);
        }
    }
#else
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
    //Iterate through the reverse map
    auto keys = _reversepossibleBindingsstencilvec[idx][idx2][t];//keys that contain t as
    // value in possiblebindings
    for(auto k:keys){
        //get the iterator range that corresponds to this key.
        auto range = _possibleBindingsstencilvec[idx][idx2].equal_range(k);
        //iterate through the range
        for(auto it = range.first; it != range.second;){
            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                _possibleBindingsstencilvec[idx][idx2].erase(it++);
            }
            else ++it;
        }
    }
    //remove from the reverse map.
    _reversepossibleBindingsstencilvec[idx][idx2][t].clear();

    int newN = _possibleBindingsstencilvec[idx][idx2].size();
    fManagervec[idx][idx2]->updateBindingReaction(newN);
    //remove all neighbors which have this binding site pair
    //Go through enclosing compartments to remove all entries with the current cylinder
    // and binding site as values.
    for(auto nc: _compartment->getenclosingNeighbours()){
        if(nc != _compartment) {
            auto m = nc->getHybridBindingSearchManager();

            //Iterate through the reverse map
            auto keys = m->_reversepossibleBindingsstencilvec[idx][idx2][t];//keys that
            // contain t as
            // value in possiblebindings
            for(auto k:keys){
                //get the iterator range that corresponds to this key.
                auto range = m->_possibleBindingsstencilvec[idx][idx2].equal_range(k);
                //iterate through the range
                for(auto it = range.first; it != range.second;){
                    if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                        m->_possibleBindingsstencilvec[idx][idx2].erase(it++);
                    }
                    else ++it;
                }
            }
            //remove from the reverse map.
            m->_reversepossibleBindingsstencilvec[idx][idx2][t].clear();

//        int oldN = m->_possibleBindingsstencilvec[idx][idx2].size();

/*            for (auto it = m->_possibleBindingsstencilvec[idx][idx2].begin(); it !=
                 m->_possibleBindingsstencilvec[idx][idx2].end();) {

                if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                    m->_possibleBindingsstencilvec[idx][idx2].erase(it++);

                else ++it;
            }*/
            int newNOther = m->_possibleBindingsstencilvec[idx][idx2].size();
            m->fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }
#endif
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
        }
    }
}

void HybridBindingSearchManager::checkoccupancySIMD(short idvec[2]){

    short idx = idvec[0];
    short idx2 = idvec[1];
    short _filamentType = 0;
    short _complimentaryfilamentType = 0;
    auto speciesboundvec = SysParams::Mechanics().speciesboundvec;
    auto maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    vector<int> CIDvec(Cylinder::vectormaxsize);
    auto ccylindervec = CUDAcommon::serlvars.ccylindervec;
    short bstatepos = bstateposvec[idx][idx2];

    for(auto cyl: Cylinder::getCylinders())
        CIDvec[cyl->_dcIndex] = cyl->getID();

    auto pbs = googlepossible[idx][idx2];

    for(auto pair = pbs.begin(); pair != pbs.end(); pair++){

        //Key
        uint32_t leg1 = pair->first;
        vector<uint32_t> leg2 = pair->second;

        uint32_t cID1 = leg1 >> 4;
        uint32_t bsite1 = mask & leg1;

        auto iter = find(CIDvec.begin(), CIDvec.end(), cID1);
        uint cIndex1 = iter - CIDvec.begin();
        CCylinder* ccyl1 = ccylindervec[cIndex1];

        short it1 = SysParams::Chemistry().bindingSites[_filamentType][bsite1];

        bool boundstate1 = false;
        bool vecboundstate1 = speciesboundvec[bstatepos][maxnbs * cIndex1 + bsite1];
        if(bstatepos == 1) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);

        }
        if(bstatepos == 2) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
        }


        //Values
        for(auto V:leg2){
            uint32_t cID2 = V >> 4;
            uint32_t bsite2 = mask & V;
            auto iter = find(CIDvec.begin(), CIDvec.end(), cID2);
            uint cIndex2 = iter - CIDvec.begin();
            CCylinder* ccyl2 = ccylindervec[cIndex2];

            short it2 = SysParams::Chemistry().bindingSites[_complimentaryfilamentType][bsite2];

            bool boundstate2 = false;
            bool vecboundstate2 = speciesboundvec[bstatepos][maxnbs * cIndex2 + bsite2];
            if(bstatepos == 1) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }
            if(bstatepos == 2) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }

            if(boundstate1 == false || boundstate2 == false) {

                std::cout << "OOPS occupied species exist " << boundstate1 << " "
                          << boundstate2 <<" "<<vecboundstate1<<" "<<vecboundstate2<< endl;
                cout<<"bstate pos 1-Linker, 2-Motor "<<bstatepos<<endl;

                std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<
                                   _compartment->coordinates()[1]<<" "<<
                                   _compartment->coordinates()[2]<<" L/M "<<
                                    ccyl1->getCylinder()->getID()<<" "<<it1<<" "<<
                                    ccyl2->getCylinder()->getID()<<" "<<it2<<" "<<cID1<<" "<<
                                    cID2<<endl;
//                exit(EXIT_FAILURE);

            }

        }
    }
}

void HybridBindingSearchManager::checkoccupancySIMD(short idvec[2], gdmap& densehashmap){

    short idx = idvec[0];
    short idx2 = idvec[1];
    short _filamentType = 0;
    short _complimentaryfilamentType = 0;
    auto speciesboundvec = SysParams::Mechanics().speciesboundvec;
    auto maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    vector<int> CIDvec(Cylinder::vectormaxsize);
    auto ccylindervec = CUDAcommon::serlvars.ccylindervec;
    short bstatepos = bstateposvec[idx][idx2];

    for(auto cyl: Cylinder::getCylinders())
        CIDvec[cyl->_dcIndex] = cyl->getID();

    auto pbs = densehashmap;

    for(auto pair = pbs.begin(); pair != pbs.end(); pair++){

        //Key
        uint32_t leg1 = pair->first;
        vector<uint32_t> leg2 = pair->second;

        uint32_t cID1 = leg1 >> 4;
        uint32_t bsite1 = mask & leg1;

        auto iter = find(CIDvec.begin(), CIDvec.end(), cID1);
        uint cIndex1 = iter - CIDvec.begin();
        CCylinder* ccyl1 = ccylindervec[cIndex1];

        short it1 = SysParams::Chemistry().bindingSites[_filamentType][bsite1];

        bool boundstate1 = false;
        bool vecboundstate1 = speciesboundvec[bstatepos][maxnbs * cIndex1 + bsite1];
        if(bstatepos == 1) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);

        }
        if(bstatepos == 2) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
        }


        //Values
        for(auto V:leg2){
            uint32_t cID2 = V >> 4;
            uint32_t bsite2 = mask & V;
            auto iter = find(CIDvec.begin(), CIDvec.end(), cID2);
            uint cIndex2 = iter - CIDvec.begin();
            CCylinder* ccyl2 = ccylindervec[cIndex2];

            short it2 = SysParams::Chemistry().bindingSites[_complimentaryfilamentType][bsite2];

            bool boundstate2 = false;
            bool vecboundstate2 = speciesboundvec[bstatepos][maxnbs * cIndex2 + bsite2];
            if(bstatepos == 1) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }
            if(bstatepos == 2) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }

            if(boundstate1 == false || boundstate2 == false) {

                std::cout << "OOPS occupied species exist " << boundstate1 << " "
                          << boundstate2 <<" "<<vecboundstate1<<" "<<vecboundstate2<< endl;
                cout<<"bstate pos 1-Linker, 2-Motor "<<bstatepos<<endl;

                std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<
                         _compartment->coordinates()[1]<<" "<<
                         _compartment->coordinates()[2]<<" L/M "<<
                         ccyl1->getCylinder()->getID()<<" "<<it1<<" "<<
                         ccyl2->getCylinder()->getID()<<" "<<it2<<" "<<cID1<<" "<<
                         cID2<<endl;
                exit(EXIT_FAILURE);

            }

        }
    }
}


void HybridBindingSearchManager::updateAllPossibleBindingsstencilHYBD() {
    for (int idx = 0; idx < totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (int idx2 = 0; idx2 < countbounds; idx2++) {
            _possibleBindingsstencilvec[idx][idx2].clear();
            _reversepossibleBindingsstencilvec[idx][idx2].clear();
        }
    }

    double min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    double minveca[2];
    double maxveca[2];
    double* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    double *coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    int Ncylincmp = _compartment->getCylinders().size();
    int cindexvec[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int> ncindex; //helper vector


    int Ncyl = Cylinder::getCylinders().size();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int idx; int idx2;
//Go through all filament types in our simulation
    for(idx = 0; idx<totaluniquefIDpairs;idx++) {

        long id = 0;
        ncindices.clear();
        auto fpairs = _filamentIDvec[idx].data();
//        int offset1 = SysParams::Mechanics().bsoffsetvec.at(fpairs[0]);
//        int offset2 = SysParams::Mechanics().bsoffsetvec.at(fpairs[1]);
        int offset1 = 0;
        int offset2 = 0;
        int nbs1 = SysParams::Chemistry().bindingSites[fpairs[0]].size();
        int nbs2 = SysParams::Chemistry().bindingSites[fpairs[1]].size();
        int Nbounds = _rMaxsqvec[idx].size();
//        auto bindingsitevec1 = SysParams::Chemistry().bindingSites[fpairs[0]];
//        auto bindingsitevec2 = SysParams::Chemistry().bindingSites[fpairs[1]];

        //Go through cylinders in the compartment, get the half of neighbors.
        for (auto c : _compartment->getCylinders()) {
            cindexvec[id] = c->_dcIndex;
            id++;
            auto Neighbors = _HneighborList->getNeighborsstencil(HNLIDvec[idx], c);
            ncindex.reserve(Neighbors.size());
            for (auto cn : Neighbors) {
                if (c->getID() > cn->getID())
                    ncindex.push_back(cn->_dcIndex);
            }
            ncindices.push_back(ncindex);
            ncindex.clear();
        }

        //Go through cylinder
        for (int i = 0; i < Ncylincmp; i++) {

            int cindex = cindexvec[i];
            short complimentaryfID;
            double x1[3], x2[3];
            double X1X2[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
            int *cnindices = ncindices[i].data();
            cylinder c = cylindervec[cindex];

            if (c.type != fpairs[0] && c.type != fpairs[1]) continue;
            else if(c.type == fpairs[0]) complimentaryfID = fpairs[1];
            else complimentaryfID = fpairs[0];

            memcpy(x1, &coord[3 * c.bindices[0]], 3 * sizeof(double));
            memcpy(x2, &coord[3 * c.bindices[1]], 3 * sizeof(double));

            //Go through the neighbors of the cylinder
            for (int arraycount = 0; arraycount < ncindices[i].size(); arraycount++) {

                int cnindex = cnindices[arraycount];
                cylinder cn = cylindervec[cnindex];

//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
//              not contain ncs that will fail this cndn.
                if (c.filamentID == cn.filamentID) continue;
                if(c.type != complimentaryfID) continue;

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

                            if (status1 && status2) continue;

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
                                    minveca[1] > maxparamcyl2[idx][1]) continue;
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
                                    maxveca[1] < minparamcyl2[idx][1]) continue;
                            }

                            for (int pos2 = 0; pos2 < nbs2; pos2++) {

                                if (areEqual(boundstate[bstatepos][offset2 + maxnbs * cn.cindex + pos2], 1.0)) {

                                    //check distances..
                                    auto mp2 = bindingsites2[idx][pos2];
                                    if (!status2)
                                        if (mp2 < maxveca[0] || mp2 > maxveca[1]) continue;
                                    if (!status1)
                                        if (mp2 > minveca[0] && mp2 < minveca[1]) continue;


                                    auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][pos1];
                                    auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][pos2];


                                    auto t1 = tuple<CCylinder *, short>(ccylvec[cindex],
                                                                        it1);
                                    auto t2 = tuple<CCylinder *, short>(ccylvec[cnindex],
                                                                        it2);

                                    //add in correct order
                                    _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                                    _reversepossibleBindingsstencilvec[idx][idx2][t2]
                                            .push_back(t1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //Place in the appropriate BindingManager
    for(short idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (short idx2 = 0; idx2 < countbounds; idx2++) {
            int newNOther = _possibleBindingsstencilvec[idx][idx2].size();
            fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }
}

void HybridBindingSearchManager::countNpairsfound(short idvec[2]){
    short idx = idvec[0];
    short idx2 = idvec[1];
    int N = 0;
    Nbindingpairs[idx][idx2] = 0;
    for (auto iter = googlepossible[idx][idx2].begin(); iter !=
             googlepossible[idx][idx2].end(); iter++) {
        N += iter->second.size();
    }
    Nbindingpairs[idx][idx2] = N;
}

void HybridBindingSearchManager::updateAllPossibleBindingsstencil() {

#ifdef SIMDBINDINGSEARCH
if(true) {
    uint dim = Totallinkermotor;
    if (dim == 1) {
        calculatebspairsself<1, true>(bspairslinkerself);
        calculatebspairsenclosed<1, false>(bspairslinker);
        calculatebspairsself<1, true>(bspairsmotorself);
        calculatebspairsenclosed<1, false>(bspairsmotor);
    } else if (dim == 2) {
        calculatebspairsself<2, true>(bspairs2self);
        calculatebspairsenclosed<2, false>(bspairs2);
    }
/*    else if(dim ==3){ template specializations not defined in dist_avx_par.cpp
        calculatebspairsself<3,true>(bspairs3self);
        calculatebspairsenclosed<3,false>(bspairs3);
    }*/
        /*else if(dim ==4){
            calculatebspairsself<4,true>();
            calculatebspairsenclosed<4,false>();
        }*/
    else {
        cout << "Number of linker/motor pairs mentioned in chemistryinput.txt exceeds "
                "maximum pairs allowed. Exiting.." << endl;
        exit(EXIT_FAILURE);
    }
/*    for(uint idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (uint idx2 = 0; idx2 < countbounds; idx2++) {
            uint newNOther = _possibleBindingsstencilSIMDvec[idx][idx2].size();
            std::cout<<"SIMD "<<newNOther<<endl;
//            fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }*/
//#else
}
#endif
#ifdef HYBRID_NLSTENCILLIST
if(false) {
    /*auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    int idx2 = 1;
    _possibleBindingsstencilvec[0][1].clear();
    int idx = 0;
    for (auto c : _compartment->getCylinders()) {
        auto fpairs = _filamentIDvec[idx].data();
        short _filamentType = fpairs[0];
        if (c->getType() != _filamentType) continue;

        auto cc = c->getCCylinder();

        for (auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
             it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {

            //now re add valid binding sites
            if (areEqual(cc->getCMonomer(*it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(),
                         1.0)) {

                for (auto cn:_HneighborList->getNeighborsstencil(HNLIDvec[idx], c)) {
                    //loop through neighbors
                    //now re add valid based on CCNL

                    if (cn->getParent() == c->getParent()) continue;
                    if (cn->getType() != _filamentType) continue;

                    auto ccn = cn->getCCylinder();

                    for (auto it2 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                         it2 !=
                         SysParams::Chemistry().bindingSites[_filamentType].end(); it2++) {

                        if (areEqual(ccn->getCMonomer(*it2)->speciesBound(
                                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(),
                                     1.0)) {

                            //check distances..
                            auto mp1 = (float) *it1 /
                                       SysParams::Geometry().cylinderNumMon[_filamentType];
                            auto mp2 = (float) *it2 /
                                       SysParams::Geometry().cylinderNumMon[_filamentType];

                            auto x1 = c->getFirstBead()->coordinate;
                            auto x2 = c->getSecondBead()->coordinate;
                            auto x3 = cn->getFirstBead()->coordinate;
                            auto x4 = cn->getSecondBead()->coordinate;

                            auto m1 = midPointCoordinate(x1, x2, mp1);
                            auto m2 = midPointCoordinate(x3, x4, mp2);

                            double dist = twoPointDistancesquared(m1, m2);

                            if (dist > _rMaxsqvec[idx][idx2] || dist <
                                                                _rMinsqvec[idx][idx2])
                                continue;

                            auto t1 = tuple<CCylinder *, short>(cc, *it1);
                            auto t2 = tuple<CCylinder *, short>(ccn, *it2);

                            //add in correct order
                            if (c->getID() > cn->getID()) {
                                _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                            }
                        }
                    }
                }
            }
        }
    }


    for (int idx = 0; idx < totaluniquefIDpairs; idx++) {
        int countbounds = _rMaxsqvec[idx].size();
        for (int idx2 = 0; idx2 < countbounds; idx2++) {
            _possibleBindingsstencilvec[idx][idx2].clear();
            _reversepossibleBindingsstencilvec[idx][idx2].clear();
        }
    }*/

}
if(false) {
    /*double min1, min2, max1, max2;
    bool status1 = true;
    bool status2 = true;
    double minveca[2];
    double maxveca[2];
    double *cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    double *coord = CUDAcommon::getSERLvars().coord;
    auto cylindervec = CUDAcommon::getSERLvars().cylindervec;
    int Ncylincmp = _compartment->getCylinders().size();
    int cindexvec[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int> ncindex; //helper vector


    int Ncyl = Cylinder::getCylinders().size();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int idx;
    int idx2;
//Go through all filament types in our simulation
    for (idx = 0; idx < totaluniquefIDpairs; idx++) {

        long id = 0;
        ncindices.clear();
        auto fpairs = _filamentIDvec[idx].data();
        int offset1 = SysParams::Mechanics().bsoffsetvec.at(fpairs[0]);
        int offset2 = SysParams::Mechanics().bsoffsetvec.at(fpairs[1]);
        int offset1 = 0;
        int offset2 = 0;
        int nbs1 = SysParams::Chemistry().bindingSites[fpairs[0]].size();
        int nbs2 = SysParams::Chemistry().bindingSites[fpairs[1]].size();
        int Nbounds = _rMaxsqvec[idx].size();
//        auto bindingsitevec1 = SysParams::Chemistry().bindingSites[fpairs[0]];
//        auto bindingsitevec2 = SysParams::Chemistry().bindingSites[fpairs[1]];

        //Go through cylinders in the compartment, get the half of neighbors.
        for (auto c : _compartment->getCylinders()) {
            cindexvec[id] = c->_dcIndex;
            id++;
            auto Neighbors = _HneighborList->getNeighborsstencil(HNLIDvec[idx], c);
            ncindex.reserve(Neighbors.size());
            for (auto cn : Neighbors) {
                if (c->getID() > cn->getID())
                    ncindex.push_back(cn->_dcIndex);
            }
            ncindices.push_back(ncindex);
            ncindex.clear();
        }

        //Go through cylinder
        for (int i = 0; i < Ncylincmp; i++) {

            int cindex = cindexvec[i];
            short complimentaryfID;
            double x1[3], x2[3];
            double X1X2[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
            int *cnindices = ncindices[i].data();
            cylinder c = cylindervec[cindex];

            if (c.type != fpairs[0] && c.type != fpairs[1]) continue;
            else if (c.type == fpairs[0]) complimentaryfID = fpairs[1];
            else complimentaryfID = fpairs[0];

            memcpy(x1, &coord[3 * c.bindices[0]], 3 * sizeof(double));
            memcpy(x2, &coord[3 * c.bindices[1]], 3 * sizeof(double));

            //Go through the neighbors of the cylinder
            for (int arraycount = 0; arraycount < ncindices[i].size(); arraycount++) {

                int cnindex = cnindices[arraycount];
                cylinder cn = cylindervec[cnindex];

//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
//              not contain ncs that will fail this cndn.
                if (c.filamentID == cn.filamentID) continue;
                if (c.type != complimentaryfID) continue;

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
                        if (areEqual(
                                boundstate[bstatepos][offset1 + maxnbs * c.cindex + pos1],
                                1.0)) {

                            auto mp1 = bindingsites1[idx][pos1];
                            double A = X3X4squared;
                            double B = 2.0 * X1X3dotX3X4 - 2.0 * mp1 * X3X4dotX1X2;
                            double C = X1X3squared + mp1 * mp1 * X1X2squared -
                                       2.0 * mp1 * X1X3dotX1X2;
                            double Bsq = B * B;
                            double C1 = C - _rMinsqvec[idx][idx2];
                            double C2 = C - _rMaxsqvec[idx][idx2];
                            double b2m4ac1 = Bsq - 4 * A * C1;
                            double b2m4ac2 = Bsq - 4 * A * C2;


                            status1 = b2m4ac1 < 0;
                            status2 = b2m4ac2 < 0;

                            if (status1 && status2) continue;

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
                                    minveca[1] > maxparamcyl2[idx][1])
                                    continue;
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
                                    maxveca[1] < minparamcyl2[idx][1])
                                    continue;
                            }

                            for (int pos2 = 0; pos2 < nbs2; pos2++) {

                                if (areEqual(
                                        boundstate[bstatepos][offset2 + maxnbs * cn.cindex +
                                                              pos2], 1.0)) {

                                    //check distances..
                                    auto mp2 = bindingsites2[idx][pos2];
                                    if (!status2)
                                        if (mp2 < maxveca[0] || mp2 > maxveca[1]) continue;
                                    if (!status1)
                                        if (mp2 > minveca[0] && mp2 < minveca[1]) continue;


                                    auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][pos1];
                                    auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][pos2];


                                    auto t1 = tuple<CCylinder *, short>(ccylvec[cindex],
                                                                        it1);
                                    auto t2 = tuple<CCylinder *, short>(ccylvec[cnindex],
                                                                        it2);

                                    //add in correct order
                                    _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                                    _reversepossibleBindingsstencilvec[idx][idx2][t2]
                                            .push_back(t1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }*/
}
#endif

    //Place in the appropriate BindingManager
/*    for(short idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (short idx2 = 0; idx2 < countbounds; idx2++) {
            int newNOther = _possibleBindingsstencilvec[idx][idx2].size();
            fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }*/
}

template <uint D, bool SELF>
void HybridBindingSearchManager::calculatebspairsself(dist::dOut<D, SELF>& bspairsoutSself){
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    auto cylcmp1 = _compartment->Cyldcindexvec;
    short bstatepos;
//    int i = 0;

    //@{ SELF
//    std::cout<<++i<<" Cmp IDs "<<_compartment->getID()<<" "<<_compartment->getID()<<endl;
//    bspairsoutSself.init_dout(10000, {900.0f, 1600.0f, 30625.0f, 50625.0f});

    minsfind = chrono::high_resolution_clock::now();
    bspairsoutSself.reset_counters();
    dist::find_distances(bspairsoutSself, _compartment->bscoords, t_avx_par);

    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runfind(minefind - minsfind);
    findtime += elapsed_runfind.count();
    //Loop through output
    short idx = 0;
    short idx2 = 0;
    short Nbounds = _rMaxsqvec[idx].size();
    for(uint dim = 0; dim < D;dim++){
        bstatepos= bstateposvec[idx][idx2];
        uint N = bspairsoutSself.counter[dim];
        for(uint pid = 0; pid < N; pid++) {
            uint16_t site1 = bspairsoutSself.dout[2 * dim][pid];
            uint16_t site2 = bspairsoutSself.dout[2 * dim + 1][pid];

            uint16_t cmpcylidx1 = site1;
            uint16_t cmpcylidx2 = site2;
            cmpcylidx1 = (cmpcylidx1 >> 4);
            cmpcylidx2 = (cmpcylidx2 >> 4);
            int cIndex1 = cylcmp1[cmpcylidx1];
            int cIndex2 = cylcmp1[cmpcylidx2];
            auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
            auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];
            short bsite1 = mask & site1;
            short bsite2 = mask & site2;

            auto fpairs = _filamentIDvec[idx].data();
            //Check if other conditions are satisfied
            if(cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
            .filamentposition - cylinder2.filamentposition) <=2)
                continue;

            if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] || cylinder1
               .type == fpairs[1] && cylinder2.type == fpairs[0]))
            {
                if(boundstate[bstatepos][maxnbs * cIndex1 + bsite1] &&
                    boundstate[bstatepos][maxnbs * cIndex2 + bsite2])
                {
                    if (cIndex1 > cIndex2) {
                        auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                        auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];
                        auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                        auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);
                        _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                        _reversepossibleBindingsstencilvec[idx][idx2][t2].push_back(t1);
                    } else {
                        auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                        auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];
                        auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                        auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);
                        _possibleBindingsstencilvec[idx][idx2].emplace(t2, t1);
                        _reversepossibleBindingsstencilvec[idx][idx2][t1].push_back(t2);
                    }
                }
            }
        }

        if(idx2 == Nbounds - 1){
            idx++;
            idx2 = 0;
            Nbounds = _rMaxsqvec[idx].size();
        } else
            idx2++;

    }
}

template<uint D, bool SELF>
void HybridBindingSearchManager::calculatebspairsenclosed (dist::dOut<D,SELF>& bspairsoutS){

    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    auto cylcmp1 = _compartment->Cyldcindexvec;
    short bstatepos;
    int i =1;
//    dist::dOut<D, SELF> bspairsoutS = bspairs<D, SELF>();
//    = bspairs<D,SELF>();

    for(auto ncmp: _compartment->getuniquepermuteNeighbours()){
        auto cylcmp2 = ncmp->Cyldcindexvec;
        minsfind = chrono::high_resolution_clock::now();
        bspairsoutS.reset_counters();
        dist::find_distances(bspairsoutS, _compartment->bscoords, ncmp->bscoords, t_avx_par);
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runfind(minefind - minsfind);
        findtime += elapsed_runfind.count();
        //Loop through output
        short idx = 0;
        short idx2 = 0;
        short Nbounds = _rMaxsqvec[idx].size();
        for(uint dim = 0; dim < D;dim++){
            bstatepos= bstateposvec[idx][idx2];
            uint N = bspairsoutS.counter[dim];
            for(uint pid = 0; pid < N; pid++) {
                uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
                uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

                uint16_t cmpcylidx1 = site1;
                uint16_t cmpcylidx2 = site2;
                cmpcylidx1 = (cmpcylidx1 >> 4);
                cmpcylidx2 = (cmpcylidx2 >> 4);
                int cIndex1 = cylcmp1[cmpcylidx1];
                int cIndex2 = cylcmp2[cmpcylidx2];
                auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
                auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];
                short bsite1 = mask & site1;
                short bsite2 = mask & site2;

                auto fpairs = _filamentIDvec[idx].data();
                //Check if other conditions are satisfied
                //Check if other conditions are satisfied
                if(cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
                   .filamentposition - cylinder2.filamentposition) <=2)
                    continue;
                if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] || cylinder1
                     .type == fpairs[1] && cylinder2.type == fpairs[0])) {
                    if((areEqual(boundstate[bstatepos][maxnbs * cIndex1 + bsite1], 1.0)) &&
                       (areEqual(boundstate[bstatepos][maxnbs * cIndex2 + bsite2], 1.0))) {
                        //for(uint idx = 0; idx < N; idx++){
                        if (cIndex1 > cIndex2) {
                            auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                            auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];
                            auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                            auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);
                            _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                            _reversepossibleBindingsstencilvec[idx][idx2][t2].push_back(t1);
                        } else {
                            auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                            auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];
                            auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                            auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);
                            ncmp->getHybridBindingSearchManager()
                                    ->_possibleBindingsstencilvec[idx][idx2].emplace(t2, t1);
                            ncmp->getHybridBindingSearchManager()
                            ->_reversepossibleBindingsstencilvec[idx][idx2][t1]
                                .push_back(t2);
                        }
                    }
                    //}
                }
            }

            if(idx2 == Nbounds - 1){
                idx++;
                idx2 = 0;
                Nbounds = _rMaxsqvec[idx].size();
            } else
                idx2++;
        }
    }
}

void HybridBindingSearchManager::updateAllPossibleBindingsstencilSIMDV2() {

#ifdef SIMDBINDINGSEARCH
    short idvecL[2] = {0,0};
    short idvecM[2] = {0,1};

    calculatebspairsLMself<1,true, true>(bspairslinkerself, idvecL);
    calculatebspairsLMenclosed<1,false, true>(bspairslinker,idvecL);

    calculatebspairsLMself<1,true, false>(bspairsmotorself, idvecM);
    calculatebspairsLMenclosed<1,false, false>(bspairsmotor,idvecM);

   /* cout<<"Size of parsed vectors "<<getpairsLinkerorMotor<true>().size()
            <<" "<<getpairsLinkerorMotor<false>().size()<<endl;*/

    //Get positions (positions based on cID_bs and nCID_bs) of keys and values in
    // pairsLinker and pairsMotor vectors
/*    vector<vector<uint>> positionreferencesL, positionreferencesM;
    vector<vector<uint32_t>> referencevector = {cID_bs, ncID_bs};
    vector<uint> temp;
    positionreferencesL.push_back(temp);
    positionreferencesL.push_back(temp);
    positionreferencesM.push_back(temp);
    positionreferencesM.push_back(temp);

    findIDsincylinderIDvector<1>(positionreferencesL, referencevector);
    findIDsincylinderIDvector<0>(positionreferencesM, referencevector);*/

    singlepassparseSIMDout<1,true>(idvecL);
    singlepassparseSIMDout<1,false>(idvecM);

    /*if(false) {
        //STEP 1 get unique keys
        vector<uint32_t> cID_bs = _compartment->cID_bs;
        vector<uint32_t> ncID_bs;

        //STEP 2 Sort keys & reverse keys
        minsfind = chrono::high_resolution_clock::now();
        thrust::stable_sort(thrust::omp::par, cID_bs.begin(), cID_bs.end());
        ncID_bs = cID_bs;
        for(auto ncmp: _compartment->getuniquepermuteNeighbours()){
            ncID_bs.insert(ncID_bs.end(),ncmp->cID_bs.begin(),ncmp->cID_bs.end());
        }
        thrust::stable_sort(thrust::omp::par, ncID_bs.begin(), ncID_bs.end());
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_sortkys(minefind - minsfind);
        std::cout << "Time sort unique keys " << elapsed_sortkys.count() << endl;

        //STEP 3 write linkerpairs and motorpairs in two vectors
        vector<vector<uint32_t>> pairlinker2D;
        vector<vector<uint32_t>> pairmotor2D;
        vector<vector<uint32_t>> Rpairlinker2D;
        vector<vector<uint32_t>> Rpairmotor2D;
        vector<uint32_t> leg1, leg2;
        //Linker
        uint N = pairslinker.size()/2;
        leg1.resize(N);
        leg2.resize(N);
        for(int i = 0; i < N; i++){
            leg1[i] = pairslinker[2*i];
            leg2[i] = pairslinker[2*i + 1];
        }
        pairlinker2D.push_back(leg1);
        pairlinker2D.push_back(leg2);
        Rpairlinker2D.push_back(leg2);
        Rpairlinker2D.push_back(leg1);
        //Motor
        N = pairsmotor.size()/2;
        leg1.resize(N);
        leg2.resize(N);
        for(int i = 0; i < N; i++){
            leg1[i] = pairsmotor[2*i];
            leg2[i] = pairsmotor[2*i + 1];
        }
        pairmotor2D.push_back(leg1);
        pairmotor2D.push_back(leg2);
        Rpairmotor2D.push_back(leg2);
        Rpairmotor2D.push_back(leg1);

        //SORT LINKER
        minsfind = chrono::high_resolution_clock::now();
        thrust::stable_sort_by_key(thrust::omp::par, pairlinker2D[0].begin(), pairlinker2D[0]
         .end(), pairlinker2D[1].begin());
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_stablethrust(minefind - minsfind);
        std::cout << "Lthrust stable sort parallel " << elapsed_stablethrust.count() << endl;

        //Lower Bound Linker
        vector<int> Lcountveckey;
        Lcountveckey.resize(cID_bs.size());

        minsfind = chrono::high_resolution_clock::now();
        thrust::lower_bound(thrust::omp::par, pairlinker2D[0].begin(), pairlinker2D[0].end(),
                    cID_bs.begin(), cID_bs.end(), Lcountveckey.begin());
        minefind = chrono::high_resolution_clock::now();

*//*        for(auto L:pairlinker2D[0])
            cout<<L<<" ";
        cout<<endl;

        for(auto L:cID_bs)
            cout<<L<<" ";
        cout<<endl;

        for(auto L:Lcountveckey)
            cout<<L<<" ";
        cout<<endl;

        cout<<Lcountveckey.size()<<endl;*//*
        chrono::duration<double> elapsed_countkey(minefind - minsfind);
        std::cout << "Lthrust countKey " << elapsed_countkey.count() << endl;

        //SORT REVERSE LINKER
        minsfind = chrono::high_resolution_clock::now();
        thrust::stable_sort_by_key(thrust::omp::par, Rpairlinker2D[0].begin(),
                                   Rpairlinker2D[0].end(), Rpairlinker2D[1].begin());
        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_stableRthrust(minefind - minsfind);
        std::cout << "LRthrust stable sort parallel " <<elapsed_stableRthrust.count() <<endl;

        //Lower Bound Reverse Linker
        vector<int> LcountvecRkey;
        LcountvecRkey.resize(ncID_bs.size());

        minsfind = chrono::high_resolution_clock::now();
        thrust::lower_bound(thrust::omp::par, Rpairlinker2D[0].begin(), Rpairlinker2D[0]
        .end(), ncID_bs.begin(), ncID_bs.end(), LcountvecRkey.begin());
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_countRkey(minefind - minsfind);
        std::cout << "LRthrust countKey " << elapsed_countRkey.count() << endl;

        //SORT MOTOR
        minsfind = chrono::high_resolution_clock::now();
        thrust::stable_sort_by_key(thrust::omp::par, pairmotor2D[0].begin(), pairmotor2D[0]
                .end(), pairmotor2D[1].begin());
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_stableMthrust(minefind - minsfind);
        std::cout << "Mthrust stable sort parallel " << elapsed_stableMthrust.count()<<endl;

        //Lower Bound Motor
        vector<int> Mcountveckey;
        Mcountveckey.resize(cID_bs.size());

        minsfind = chrono::high_resolution_clock::now();
        thrust::lower_bound(thrust::omp::par, pairmotor2D[0].begin(), pairmotor2D[0].end(),
                            cID_bs.begin(), cID_bs.end(), Mcountveckey.begin());
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_countMkey(minefind - minsfind);
        std::cout << "Mthrust countKey " << elapsed_countMkey.count() << endl;

        //SORT REVERSE MOTOR
        minsfind = chrono::high_resolution_clock::now();
        thrust::stable_sort_by_key(thrust::omp::par, Rpairmotor2D[0].begin(),
                Rpairmotor2D[0].end(), Rpairmotor2D[1].begin());
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_stableMRthrust(minefind - minsfind);
        std::cout << "MRthrust stable sort parallel " << elapsed_stableMRthrust.count()
        <<endl;

        //Lower Bound REVERSE Motor
        vector<int> McountvecRkey;
        McountvecRkey.resize(ncID_bs.size());

        minsfind = chrono::high_resolution_clock::now();
        thrust::lower_bound(thrust::omp::par, Rpairmotor2D[0].begin(), Rpairmotor2D[0]
        .end(), ncID_bs.begin(), ncID_bs.end(), McountvecRkey.begin());
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_countMRkey(minefind - minsfind);
        std::cout << "MRthrust countKey " << elapsed_countMRkey.count() << endl;

        dense_hash_map<uint32_t , vector<uint>, hash<uint32_t>> Lmap, LRmap, Mmap, MRmap;
        Lmap.set_empty_key(4294967295);
        Lmap.set_deleted_key(4294967294);
        LRmap.set_empty_key(4294967295);
        LRmap.set_deleted_key(4294967294);
        Mmap.set_empty_key(4294967295);
        Mmap.set_deleted_key(4294967294);
        MRmap.set_empty_key(4294967295);
        MRmap.set_deleted_key(4294967294);
        uint Lfirst, Lcount, Llast, Mfirst, Mcount, Mlast;

        minsfind = chrono::high_resolution_clock::now();
        for(uint i = 0; i < cID_bs.size(); i++){
            uint32_t key = cID_bs[i];
            //Linker
            Lfirst = Lcountveckey[i];
            Llast = Lcountveckey[i+1];
            Lcount = Llast - Lfirst;
            vector<uint> value = {Lfirst,Llast, Lcount};
            Lmap[key].insert(Lmap[key].end(),value.begin(),value.end());
            //Motor
            Mfirst = Mcountveckey[i];
            Mlast = Mcountveckey[i+1];
            Mcount = Mlast - Mfirst;
            vector<uint> Mvalue = {Mfirst,Mlast, Mcount};
            Mmap[key].insert(Mmap[key].end(),Mvalue.begin(),Mvalue.end());
        }

        for(uint i = 0; i < ncID_bs.size(); i++){
            uint32_t key = cID_bs[i];
            //Linker
            Lfirst = LcountvecRkey[i];
            Llast = LcountvecRkey[i+1];
            Lcount = Llast - Lfirst;
            vector<uint> value = {Lfirst,Llast, Lcount};
            LRmap[key].insert(LRmap[key].end(),value.begin(),value.end());
            //Motor
            Mfirst = McountvecRkey[i];
            Mlast = McountvecRkey[i+1];
            Mcount = Mlast - Mfirst;
            vector<uint> Mvalue = {Mfirst,Mlast, Mcount};
            MRmap[key].insert(MRmap[key].end(),Mvalue.begin(),Mvalue.end());
        }
        minefind = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed_map(minefind - minsfind);
        std::cout << "maptime " << elapsed_map.count() << endl;
    }*/

    //Place in the appropriate BindingManager
/*    for(short idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (short idx2 = 0; idx2 < countbounds; idx2++) {
            int newNOther = _possibleBindingsstencilvec[idx][idx2].size();
            fManagervec[idx][idx2]->updateBindingReaction(newNOther);
        }
    }*/
//Calculate N
    minsfind = chrono::high_resolution_clock::now();
    for(short idx = 0; idx<totaluniquefIDpairs; idx++) {
        int countbounds = _rMaxsqvec[idx].size();
        for (short idx2 = 0; idx2 < countbounds; idx2++) {
            short idvec[2] = {idx, idx2};;
            countNpairsfound(idvec);
//            cout<<idx<<" "<<idx2<<" "<<Nbindingpairs[idx][idx2]<<endl;
            fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);
        }
    }
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_countsites(minefind - minsfind);
    SIMDcountbs += elapsed_countsites.count();
//    cout<<"Time to count dense hash map size "<<elapsed_countsites.count()<<endl;
#endif
}

template<bool LinkerorMotor>
void HybridBindingSearchManager::findIDsincylinderIDvector(vector<vector<uint>>& outputvector,
        vector<vector<uint32_t>>& cID_bs){

    /*minsfind = chrono::high_resolution_clock::now();
    auto searchvector = getpairsLinkerorMotor<LinkerorMotor>();
    //ID = 0 Leg1, ID = 1 Leg2
        outputvector.resize(searchvector.size());
        thrust::lower_bound(thrust::omp::par,
                            cID_bs[ID].begin(), cID_bs[ID].end(),
                            searchvector[ID].begin(), searchvector[ID].end(),
                            outputvector[ID].begin());
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_lowkys(minefind - minsfind);
    std::cout << "Time trace lower bound keys " << elapsed_lowkys.count() << endl;*/

}

void HybridBindingSearchManager::countKeys(vector<uint16_t>& keyvec, vector<unsigned int
        long>& bounds, vector<int>& countvec, vector<uint32_t>& bindingskey, uint prev,
        uint next){

    int firstkeypos = prev;
    int lastkeypos= next;
    int totalkeys = bounds[2];
    int bindingskeysize = bindingskey.size();
    int scanstartpos = 0;
    int _key = keyvec[firstkeypos];
    auto lower = std::lower_bound(bindingskey.begin(), bindingskey.end(), _key);
    scanstartpos =  lower-bindingskey.begin();
    for(int pos = firstkeypos; pos < lastkeypos; pos++){
        uint32_t _key = keyvec[pos];
        int count = 0;
        while(scanstartpos < bindingskeysize){
            if(bindingskey[scanstartpos] == _key) {
                scanstartpos++;
                count++;
            }
            else{
                countvec[pos] = count;
                break;
            }
        }
    }
}

template <uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::calculatebspairsLMself(dist::dOut<D, SELF>&
        bspairsoutSself, short idvec[2]){

    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    auto cylcmp1 = _compartment->Cyldcindexvec;

    minsfind = chrono::high_resolution_clock::now();

    bspairsoutSself.reset_counters();
    dist::find_distances(bspairsoutSself, _compartment->getSIMDcoords<LinkerorMotor>(),
            t_avx_par);

    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runfind(minefind - minsfind);
    findtimeV2 += elapsed_runfind.count();

    //MERGE INTO single vector
    //@{
    if(true) {
//        minsfind = chrono::high_resolution_clock::now();

        uint N = bspairsoutSself.counter[D-1];
/*        cout<<"Nbs SIMD "<<N<<" prev size "<<getfilID_fpospairs<LinkerorMotor>().size()
        <<" ";*/
        uint prev_size = getfilID_fpospairs<LinkerorMotor>().size();
        getpairsLinkerorMotor<LinkerorMotor>().resize
            (getpairsLinkerorMotor<LinkerorMotor>().size() + 2*N);
        getfilID_fpospairs<LinkerorMotor>().resize(getfilID_fpospairs<LinkerorMotor>()
                .size() + N);
//        cout<<"aftr size "<<getfilID_fpospairs<LinkerorMotor>().size()<<endl;
/*            getboundspeciesLinkerorMotor<LinkerorMotor>()[i].resize
            (getboundspeciesLinkerorMotor<LinkerorMotor>()[i].size() + N);*/

        //Need to append as Cylinder ID| bs and not _cIndex|bs;
       /* getpairsLinkerorMotor<LinkerorMotor>()[i].insert(
                getpairsLinkerorMotor<LinkerorMotor>()[i].end(),
                bspairsoutSself.dout[2 * dim + i].begin(),
                bspairsoutSself.dout[2 * dim + i].begin()+ N);*/
//        cout<<getpairsLinkerorMotor<LinkerorMotor>()[i].size()<<" "<<bspairsoutSself.counter[dim]<<endl;


        minsfind = chrono::high_resolution_clock::now();
        std::vector<std::thread> threads_avx;
        uint nt = nthreads;
        threads_avx.reserve(nt);
        uint prev = 0;
        uint frac = N / nt;
        uint next = frac + N % nt;
        for (uint i = 0; i < nt; ++i) {
            threads_avx.push_back(std::thread(
                    &HybridBindingSearchManager::gatherCylinderIDfromcIndex<D, SELF,
                            LinkerorMotor>, this, std::ref(bspairsoutSself), prev,
                    next, prev_size, _compartment));
            prev = next;
            next = min(N, prev + frac);
        }

        //Join
        for (auto &t : threads_avx)
            t.join();
        threads_avx.clear();

/*        for (int i = 0; i < 2; i++)
            cout << "Size " << getpairsLinkerorMotor<LinkerorMotor>()[i].size() << " " << N
                 << endl;*/

        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_append(minefind - minsfind);
//        cout << "Time parallel lookup " << elapsed_append.count() << endl;
        appendtime += elapsed_append.count();
    }
    //@}

if(false)
    parseSIMDout<D, SELF, LinkerorMotor>(bspairsoutSself, idvec, _compartment);


    //Loop through output
    if(false) {
        short idx = idvec[0];
        short idx2 = idvec[1];
        short Nbounds = _rMaxsqvec[idx].size();
        uint dim = 0;

        {
            uint N = bspairsoutSself.counter[dim];

            cout << "Size of google possible " << googlepossiblel.size() << endl;

            for (const auto &iter:googlepossiblel) {
                cout << iter.second.size() << " ";
            }
            cout << endl;
//        googlepossiblel.clear();
            //Reserve
            /*getLinkerorMotor<LinkerorMotor>(1).reserve(getLinkerorMotor<LinkerorMotor>(1).size() + N);
            getLinkerorMotor<LinkerorMotor>(2).reserve(getLinkerorMotor<LinkerorMotor>(2).size() + N);
            uint pos = getLinkerorMotor<LinkerorMotor>(1).size();*/
            for (uint pid = 0; pid < N; pid++) {
                uint16_t site1 = bspairsoutSself.dout[2 * dim][pid];
                uint16_t site2 = bspairsoutSself.dout[2 * dim + 1][pid];

                uint16_t cmpcylidx1 = site1;
                uint16_t cmpcylidx2 = site2;
                cmpcylidx1 = (cmpcylidx1 >> 4);
                cmpcylidx2 = (cmpcylidx2 >> 4);
                int cIndex1 = cylcmp1[cmpcylidx1];
                int cIndex2 = cylcmp1[cmpcylidx2];
                auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
                auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];
                short bsite1 = mask & site1;
                short bsite2 = mask & site2;

                auto fpairs = _filamentIDvec[idx].data();
                //Check if other conditions are satisfied
                if (cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
                                                                                .filamentposition -
                                                                        cylinder2.filamentposition) <=
                                                                    2)
                    continue;

/*                if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] || cylinder1
                                                                                           .type ==
                                                                                   fpairs[1] &&
                                                                                   cylinder2.type ==
                                                                                   fpairs[0])) */
                {

                    if (cIndex1 > cIndex2) {
                        /*auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                        auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];*/
                        /*auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                        auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);*/
/*                        _possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                        _reversepossibleBindingsstencilvec[idx][idx2][t2].push_back(t1);*/
                        uint32_t t1 = cIndex1 << 4 | bsite1;
                        uint32_t t2 = cIndex2 << 4 | bsite2;
//                        i++;
                        googlepossibleORIG[t1].push_back(t2);
                        googlereversepossibleORIG[t2].push_back(t1);
                        /* _possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);
                         _reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);*/

                        /*getLinkerorMotor<LinkerorMotor>(1)[Npairs] = t1;
                        getLinkerorMotor<LinkerorMotor>(2)[Npairs] = t2;
                        Npairs++;*/
                    } else {
                        /*auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                        auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];*/
                        /*auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                        auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);*/
                        /*_possibleBindingsstencilvec[idx][idx2].emplace(t2, t1);
                        _reversepossibleBindingsstencilvec[idx][idx2][t1].push_back(t2);*/
                        uint32_t t1 = cIndex1 << 4 | bsite1;
                        uint32_t t2 = cIndex2 << 4 | bsite2;
//                        i++;
                        googlepossibleORIG[t2].push_back(t1);
                        googlereversepossibleORIG[t1].push_back(t2);

/*                        _possibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);
                        _reversepossibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);*/

/*                        getLinkerorMotor<LinkerorMotor>(1)[Npairs] = t2;
                        getLinkerorMotor<LinkerorMotor>(2)[Npairs] = t1;
                        Npairs++;*/

                    }
                }
            }
        }

        cout << "Size of google possible ORIG " << googlepossibleORIG.size() << endl;
        for (const auto &iter:googlepossibleORIG) {
            cout << iter.second.size() << " ";
        }
        cout << endl;
    }
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::calculatebspairsLMenclosed (dist::dOut<D,SELF>&
        bspairsoutS, short idvec[2]){
    uint i =0;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    CCylinder **ccylvec = CUDAcommon::getSERLvars().ccylindervec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    auto cylcmp1 = _compartment->Cyldcindexvec;

    for(auto ncmp: _compartment->getuniquepermuteNeighbours()){

        minsfind = chrono::high_resolution_clock::now();

        bspairsoutS.reset_counters();
        dist::find_distances(bspairsoutS, _compartment->getSIMDcoords<LinkerorMotor>(),
                ncmp->getSIMDcoords<LinkerorMotor>(), t_avx_par);

        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runfind(minefind - minsfind);
        findtimeV2 += elapsed_runfind.count();

    //MERGE INTO single vector
    //@{
    if(true) {
        minsfind = chrono::high_resolution_clock::now();
        short dim = 0;
        uint N = bspairsoutS.counter[dim];
        uint prev_size = getfilID_fpospairs<LinkerorMotor>().size();
/*        cout<<"Nbs SIMD "<<N<<"prev size "<<getfilID_fpospairs<LinkerorMotor>().size()
            <<" ";*/
        getpairsLinkerorMotor<LinkerorMotor>().resize
                    (getpairsLinkerorMotor<LinkerorMotor>().size() + 2*N);
        getfilID_fpospairs<LinkerorMotor>().resize(getfilID_fpospairs<LinkerorMotor>()
                                                           .size() + N);
//        cout<<"aftr size "<<getfilID_fpospairs<LinkerorMotor>().size()<<endl;

            /*getboundspeciesLinkerorMotor<LinkerorMotor>()[i].resize
                    (getboundspeciesLinkerorMotor<LinkerorMotor>()[i].size() + N);*/
            /*getpairsLinkerorMotor<LinkerorMotor>()[i].insert(
                    getpairsLinkerorMotor<LinkerorMotor>()[i].end(),
                    bspairsoutS.dout[2 * dim + i].begin(),
                    bspairsoutS.dout[2 * dim + i].begin()+ N);*/
            /*cout<<getpairsLinkerorMotor<LinkerorMotor>()[i].size()<<" "<<bspairsoutS
            .counter[dim]<<endl;*/

        std::vector<std::thread> threads_avx;
        uint nt = nthreads;
        threads_avx.reserve(nt);
        uint prev = 0;
        uint frac = N / nt;
        uint next = frac + N % nt;
        for (uint i = 0; i < nt; ++i) {
            threads_avx.push_back(std::thread
                                          (&HybridBindingSearchManager::gatherCylinderIDfromcIndex<D, SELF,
                                                   LinkerorMotor>, this,
                                           std::ref(bspairsoutS), prev,
                                           next, prev_size, ncmp));
            prev = next;
            next = min(N, prev + frac);
        }

        //Join
        for (auto &t : threads_avx)
            t.join();
        threads_avx.clear();

        minefind = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_append(minefind - minsfind);
        appendtime += elapsed_append.count();
    }
    //@}

        if(false)
            parseSIMDout<D, SELF, LinkerorMotor>(bspairsoutS, idvec, ncmp);
        //Loop through output
        if(false) {
                    auto cylcmp2 = ncmp->Cyldcindexvec;
            short idx = idvec[0];
            short idx2 = idvec[1];
            short Nbounds = _rMaxsqvec[idx].size();
            uint dim = 0;
            {
                uint N = bspairsoutS.counter[dim];


                cout << "Size of google possible " << googlepossiblel.size() << endl;

                for (const auto &iter:googlepossiblel) {
                    cout << iter.second.size() << " ";
                }
                cout << endl;
//            googlepossiblel.clear();

                //reserve
                /*getLinkerorMotor<LinkerorMotor>(1).reserve(getLinkerorMotor<LinkerorMotor>(1)
                        .size() + N/2);
                getLinkerorMotor<LinkerorMotor>(2).reserve(getLinkerorMotor<LinkerorMotor>(2)
                        .size() + N/2);
                //Neighbor
                ncmp->getHybridBindingSearchManager()->getLinkerorMotor<LinkerorMotor>(1).reserve
                (ncmp->getHybridBindingSearchManager()->getLinkerorMotor<LinkerorMotor>(1).size() + N/2);
                ncmp->getHybridBindingSearchManager()->getLinkerorMotor<LinkerorMotor>(2).reserve
                (ncmp->getHybridBindingSearchManager()->getLinkerorMotor<LinkerorMotor>(2).size() + N/2);*/
//            std::cout<<"Found pairs SIMD "<<N<<" serial "<<bsparirsoutserial.counter[dim]<<endl;

                for (uint pid = 0; pid < N; pid++) {
                    uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
                    uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

                    uint16_t cmpcylidx1 = site1;
                    uint16_t cmpcylidx2 = site2;
                    cmpcylidx1 = (cmpcylidx1 >> 4);
                    cmpcylidx2 = (cmpcylidx2 >> 4);
                    int cIndex1 = cylcmp1[cmpcylidx1];
                    int cIndex2 = cylcmp2[cmpcylidx2];
                    auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
                    auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];
                    short bsite1 = mask & site1;
                    short bsite2 = mask & site2;

                    auto fpairs = _filamentIDvec[idx].data();
                    //Check if other conditions are satisfied
                    if (cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
                                                                                    .filamentposition -
                                                                            cylinder2.filamentposition) <=
                                                                        2)
                        continue;
                    if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] ||
                         cylinder1
                                 .type == fpairs[1] && cylinder2.type == fpairs[0])) {

                        //for(uint idx = 0; idx < N; idx++){
                        if (cIndex1 > cIndex2) {
/*                            auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                            auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];*/
/*                            auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                            auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);*/
                            /*_possibleBindingsstencilvec[idx][idx2].emplace(t1, t2);
                            _reversepossibleBindingsstencilvec[idx][idx2][t2]
                                    .push_back(t1);*/
                            uint32_t t1 = cIndex1 << 4 | bsite1;
                            uint32_t t2 = cIndex2 << 4 | bsite2;
                            googlepossibleORIG[t1].push_back(t2);
                            googlereversepossibleORIG[t2].push_back(t1);
//                            i++;
                            /*getLinkerorMotor<LinkerorMotor>(1)[Npairs] = (t1);
                            getLinkerorMotor<LinkerorMotor>(2)[Npairs] = (t2);
                            Npairs++;*/
/*                            _possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);
                            _reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);*/
                        } else {
/*                            auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsite1];
                            auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsite2];*/
/*                            auto t1 = tuple<CCylinder *, short>(ccylvec[cIndex1], it1);
                            auto t2 = tuple<CCylinder *, short>(ccylvec[cIndex2], it2);*/
                            /*ncmp->getHybridBindingSearchManager()
                                    ->_possibleBindingsstencilvec[idx][idx2].emplace(t2,
                                                                                     t1);
                            ncmp->getHybridBindingSearchManager()
                                    ->_reversepossibleBindingsstencilvec[idx][idx2][t1]
                                    .push_back(t2);*/
                            uint32_t t1 = cIndex1 << 4 | bsite1;
                            uint32_t t2 = cIndex2 << 4 | bsite2;
//                            i++;
                            ncmp->getHybridBindingSearchManager()
                                    ->googlepossibleORIG[t2].push_back(t1);
                            ncmp->getHybridBindingSearchManager()
                                    ->googlereversepossibleORIG[t1].push_back(t2);

/*                            uint nNpairs = ncmp->getHybridBindingSearchManager()->Npairs;
                            ncmp->getHybridBindingSearchManager()
                                    ->getLinkerorMotor<LinkerorMotor>(1)[nNpairs] = (t2);
                            ncmp->getHybridBindingSearchManager()
                                    ->getLinkerorMotor<LinkerorMotor>(2)[nNpairs] = (t1);
                            ncmp->getHybridBindingSearchManager()->Npairs++;*/

/*                            ncmp->getHybridBindingSearchManager()
                                    ->_possibleBindingsstencilvecuint[idx][idx2][t2]
                                    .push_back(t1);
                            ncmp->getHybridBindingSearchManager()
                                    ->_reversepossibleBindingsstencilvecuint[idx][idx2][t1]
                                    .emplace_back(t2);*/
                        }
                    }
                }
            }
        }
    }

}

template <uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::gatherCylinderIDfromcIndex(dist::dOut<D,SELF>&
        bspairsoutS, int first, int last, uint prev_size, Compartment* nCmp){
//    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    //Need to rewrite to give CylinderID instead.
/*    auto cylcmp1 = _compartment->Cyldcindexvec;
    auto cylcmp2 = nCmp->Cyldcindexvec;*/

    auto cylcmpcIndex1 = _compartment->Cyldcindexvec;
    auto cylcmpcIndex2 = nCmp->Cyldcindexvec;

    auto cylcmp1 = _compartment->CylcIDvec;
    auto cylcmp2 = nCmp->CylcIDvec;

    auto cylindervec = CUDAcommon::serlvars.cylindervec;

    auto cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    auto ccylindervec = CUDAcommon::serlvars.ccylindervec;

    int writeposition = 0;
    for(uint pid = first; pid < last; pid++){
        uint16_t site1 = bspairsoutS.dout[2 * (D -1)][pid];
        uint16_t site2 = bspairsoutS.dout[2 * (D - 1) + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;
        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        uint32_t cIndex1 =  cylcmpcIndex1[cmpcylidx1];
        uint32_t cIndex2 =  cylcmpcIndex2[cmpcylidx2];

        cylinder cylinder1 = cylindervec[cIndex1];
        cylinder cylinder2 = cylindervec[cIndex2];

        if(cylinder1.filamentID == cylinder2.filamentID &&
            abs(cylinder1.filamentposition - cylinder2.filamentposition) <=2){
            /*if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] ||
                 cylinder1.type == fpairs[1] && cylinder2.type == fpairs[0]))*/
            getfilID_fpospairs<LinkerorMotor>()[prev_size + pid] = 0;

        } else
            getfilID_fpospairs<LinkerorMotor>()[prev_size + pid] = 1;


        uint32_t cID1 = cylcmp1[cmpcylidx1];
        uint32_t cID2 = cylcmp2[cmpcylidx2];

        int Ncyl = Cylinder::getCylinders().size();
        if(cID1 >= Ncyl || cID2 >= Ncyl){
            cout<<"OOPS! Compartment ID does not exist "<<cID1<<" "<<cID2<<" MaxNumCyl "
                                                                           ""<<Ncyl<<endl;
        }


/*        if(cylinderpointervec[cIndex1]->getID() != cID1 ||
           cylinderpointervec[cIndex2]->getID() != cID2){
            cout<<"Found the wrong cylinder "<<cylinderpointervec[cIndex1]->getID()
                <<" "<<cID1<<" "<<cylinderpointervec[cIndex2]->getID()<<" "<<cID2<<endl;
        }

        if(cID1 != cylinder1.ID || cID2 != cylinder2.ID){
            cout<<"OOPS! Big mistake Cylinder IDs don't match "<<cID1<<" "<<cID2<<" "
            <<cylinder1.ID<<" "<<cylinder2.ID<<endl;
        }*/

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

/*        if(!LinkerorMotor) {
        short _filamentType = 0;
        short it1 = SysParams::Chemistry().bindingSites[_filamentType][bsite1];
        short it2 = SysParams::Chemistry().bindingSites[_filamentType][bsite2];
        auto ccyl1 = ccylindervec[cIndex1];
        auto ccyl2 = ccylindervec[cIndex2];
            bool boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
            bool boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
            if(boundstate1 == false || boundstate2 == false){
                cout<<"OOPS Trying to add an occupied site "<<cylinder1.ID<<" "
                    <<cylinder2.ID<<" "<<it1<<" "<<it2<<endl;
            }
        }*/

        uint32_t t1 = cID1 << 4 | bsite1;
        uint32_t t2 = cID2 << 4 | bsite2;
/*        if(bsite1>3 || bsite2>3) {
            cout << "cID " << bitset<32>(cID1) << " " << bitset<32>(cID2) << endl;
            cout << "site " << bitset<32>(bsite1) << " " << bitset<32>(bsite2) << endl;
            cout << "pack " << bitset<32>(t1) << " " << bitset<32>(t2) << endl;
            cout<<endl;
        }*/

        getpairsLinkerorMotor<LinkerorMotor>()[2 * prev_size + 2 * pid] = t1;
        getpairsLinkerorMotor<LinkerorMotor>()[2 * prev_size + 2 * pid + 1] = t2;

       /* getboundspeciesLinkerorMotor<LinkerorMotor>()[writeposition][pid] =
                boundstate[LinkerorMotor + 1][maxnbs * cIndex1 + bsite1];
        getboundspeciesLinkerorMotor<LinkerorMotor>()[writeposition + 1][pid] =
                boundstate[LinkerorMotor + 1][maxnbs * cIndex2 + bsite2];*/
    }
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::parseSIMDout
    (dist::dOut<D,SELF>& bspairsoutS, short idvec[2], Compartment* ncmp) {
    minsfind = chrono::high_resolution_clock::now();
    uint dim = 0;
    uint N = bspairsoutS.counter[dim];

    std::vector<std::thread> threads_parse;
    uint prev = 0;
    uint nt = nthreads;
    uint frac = N / nt;
    uint next = frac + N % nt;

    for (uint i = 0; i < nt; ++i) {

        if (SELF)
            threads_parse.push_back(std::thread
                                          (&HybridBindingSearchManager::threadwritepairsselfV3<D, SELF, LinkerorMotor>,
                                           this, prev, next, std::ref(bspairsoutS), idvec, i));
        else
            threads_parse.push_back(std::thread(
                    &HybridBindingSearchManager::threadwritepairsV3<D, SELF, LinkerorMotor>,
                    this, ncmp, std::ref(bspairsoutS), prev, next, idvec, i));

        prev = next;
        next = min(N, prev + frac);

    }

    //Join
    for (auto &t : threads_parse)
        t.join();
    std::vector<std::thread> threads_SIMDV2;
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD(minefind - minsfind);
    SIMDparse1 += elapsed_runSIMD.count();

    //COMMENTED OUT
    if(false){
/*        std::vector<std::thread> threads_avx;
    for (uint i = 0; i < nt; ++i) {

        if (SELF)
            threads_avx.push_back(std::thread
                                          (&HybridBindingSearchManager::threadwritepairsself<D, SELF, LinkerorMotor>,
                                           this, prev, next, std::ref(bspairsoutS), idvec,
                                           (vecmap), i));
        else
            threads_avx.push_back(std::thread(
                    &HybridBindingSearchManager::threadwritepairs<D, SELF, LinkerorMotor>,
                    this, ncmp, std::ref(bspairsoutS), vecmap, vecNmap,
                    prev, next, idvec, i));

*//*        threadwritepairs<D, SELF, LinkerorMotor>(ncmp, std::ref(bspairsoutS), vecmap,
                vecNmap, prev, next, idvec, i);*//*
*//*        threadwritepairs<D,SELF, LinkerorMotor>(prev, next, std::ref(bspairsoutS), idvec,
                                                (vecmap[i]), ncmp, i);*//*
        prev = next;
        next = min(N, prev + frac);

    }

    //Join
    for (auto &t : threads_avx)
        t.join();
    threads_avx.clear();

    minsfind = chrono::high_resolution_clock::now();
    prev = 0;
    next = frac + N % nt;
    for (uint i = 0; i < nt; ++i) {
        cout << prev << " " << next << endl;
        threadwritepairsselfV2<D, SELF, LinkerorMotor>(prev, next, std::ref(bspairsoutS),
                                                       idvec,
                                                       (vecmapV2), valuematrixvec, i);
//        threads_avx.push_back(std::thread(&HybridBindingSearchManager::threadtest,this,i));
        *//*threads_SIMDV2.push_back(std::thread
                (&HybridBindingSearchManager::threadwritepairsselfV2<D, SELF, LinkerorMotor>,
                this, prev, next, std::ref(bspairsoutS), idvec, (vecmapV2), valuematrixvec, i));*//*
        prev = next;
        next = min(N, prev + frac);
    }

    for (auto &t : threads_SIMDV2)
        t.join();
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMDV2(minefind - minsfind);
    cout << "SIMDV2 time " << elapsed_runSIMDV2.count() << " SIMDV1 time "
         << elapsed_runSIMD
                 .count() << endl;
    for (uint i = 0; i < nt; ++i) {
        cout << vecmapV2[2 * i].size() << " " << vecmapV2[2 * i + 1].size() << endl;
        cout << valuematrixvec[2 * i].size() << " " << valuematrixvec[2 * i + 1].size()
             << endl;
    }*/
    } //FALSE

    minsfind = chrono::high_resolution_clock::now();
    for(uint i=0; i<nthreads; ++i) {

        if (getvecmapLinkerorMotor<LinkerorMotor>()[2 * i].size()) {
            for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[2 * i]) {

                vector<uint32_t> sitestoappend = iter.second;
                auto mapiter = googlepossiblel.find(iter.first);
                if (mapiter == googlepossiblel.end()) {
                    sitestoappend.reserve(sitestoappend.size()+1024);
                    googlepossiblel.insert(make_pair(iter.first, sitestoappend));
                } else {
                    mapiter->second.insert(mapiter->second.end(), sitestoappend.begin(),
                            sitestoappend.end());
                }
            }
        }
        //reverse map
        if (getvecmapLinkerorMotor<LinkerorMotor>()[2 * i + 1].size()) {
            for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[2 * i + 1]) {

                vector<uint32_t> sitestoappend = iter.second;
                auto mapiter = googlereversepossiblel.find(iter.first);
                if (mapiter == googlereversepossiblel.end()) {
                    sitestoappend.reserve(sitestoappend.size()+1024);
                    googlereversepossiblel.insert(make_pair(iter.first, sitestoappend));
                } else {
                    mapiter->second.insert(mapiter->second.end(), sitestoappend.begin(),
                                           sitestoappend.end());
                }
            }
        }

        if(false) {
            /*if (!SELF) {
                if (vecNmap[2 * i].size()) {
                    for (const auto &iter : vecNmap[2 * i]) {

                        *//*vector<uint32_t> mapiter = googlepossiblel[iter.first];
                        cout<<mapiter.size()<<endl;*//*
                        vector<uint32_t> sitestoappend = iter.second;
                        auto mapiter = ncmp->getHybridBindingSearchManager()->googlepossiblel.find(
                                iter.first);
                        if (mapiter ==
                            ncmp->getHybridBindingSearchManager()->googlepossiblel.end()) {
                            sitestoappend.reserve(sitestoappend.size() + 1024);
                            ncmp->getHybridBindingSearchManager()->googlepossiblel.insert(
                                    make_pair
                                            (iter.first, sitestoappend));
                        } else {
                            mapiter->second.insert(mapiter->second.end(),
                                                   sitestoappend.begin(),
                                                   sitestoappend.end());
                        }

                    }
                }
                //reverse map
                if (vecNmap[2 * i + 1].size()) {
                    for (const auto &iter : vecNmap[2 * i + 1]) {

                        *//*vector<uint32_t> mapiter = googlepossiblel[iter.first];
                        cout<<mapiter.size()<<endl;*//*
                        vector<uint32_t> sitestoappend = iter.second;
                        auto mapiter = ncmp->getHybridBindingSearchManager()
                                ->googlereversepossiblel.find(
                                        iter.first);
                        if (mapiter ==
                            ncmp->getHybridBindingSearchManager()->googlereversepossiblel.end
                                    ()) {
                            sitestoappend.reserve(sitestoappend.size() + 1024);
                            ncmp->getHybridBindingSearchManager()->googlereversepossiblel
                                    .insert(
                                            make_pair(iter.first, sitestoappend));
                        } else {
                            mapiter->second.insert(mapiter->second.end(),
                                                   sitestoappend.begin(),
                                                   sitestoappend.end());
                        }

                    }
                }
            }*/
        } // FALSE
    }
    //Clear
    for(uint i=0; i<nthreads; ++i)
        getvecmapLinkerorMotor<LinkerorMotor>()[i].clear();
/*    if(SELF == false){
        for(uint i=0; i<nthreads; ++i)
            vecNmap[i].clear();
    }*/
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD2(minefind - minsfind);
    SIMDparse2 += elapsed_runSIMD2.count();

    //Check if there are entries
/*
    for(auto key: googlepossiblel){
        cout<<key.first<<" "<<key.second.size()<<endl;
    }
*/

    }

template<uint D, bool LinkerorMotor>
void HybridBindingSearchManager::singlepassparseSIMDout(short idvec[2]){
    minsfind = chrono::high_resolution_clock::now();
    uint dim = 0;
    uint N = getpairsLinkerorMotor<LinkerorMotor>().size()/2;

    std::vector<std::thread> threads_parse;
    uint prev = 0;
    uint nt = nthreads;
    uint frac = N / nt;
    uint next = frac + N % nt;

/*    cout<<"b4 map write checking idvec "<<idvec[0]<<" "<<idvec[1]<<endl;
    for(int i= 0; i <2*nthreads; i=i+2) {
        checkoccupancySIMD(idvec, getvecmapLinkerorMotor<LinkerorMotor>()[i]);
    }
    checkoccupancySIMD(idvec);*/

    for (uint i = 0; i < nt; ++i) {

        threads_parse.push_back(std::thread(
                    &HybridBindingSearchManager::threadedsingleparse<D, LinkerorMotor>,
                    this, prev, next, i));

        prev = next;
        next = min(N, prev + frac);

    }

    //Join
    for (auto &t : threads_parse)
        t.join();

    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD(minefind - minsfind);
    SIMDparse1 += elapsed_runSIMD.count();

    //Merge maps
    minsfind = chrono::high_resolution_clock::now();
    short Nt = nthreads;

    while(Nt>=2){
        std::vector<std::thread> threads_merge;
        threads_merge.reserve(Nt);
        for (short i = 0; i < Nt; ++i) {
            threads_merge.push_back(std::thread(
                    &HybridBindingSearchManager::mergemaps<LinkerorMotor>, this, i,Nt));
        }
        for (auto &t : threads_merge)
            t.join();
        Nt=Nt/2;
    }
    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD2(minefind - minsfind);
    SIMDparse2 += elapsed_runSIMD2.count();

    minsfind = chrono::high_resolution_clock::now();

    googlepossible[idvec[0]][idvec[1]] = getvecmapLinkerorMotor<LinkerorMotor>()[0];
    googlereversepossible[idvec[0]][idvec[1]] = getvecmapLinkerorMotor<LinkerorMotor>()[1];

    minefind = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runSIMD3(minefind - minsfind);
    SIMDparse3 += elapsed_runSIMD3.count();

    //Clear
    for(uint i=0; i<2 * nthreads; ++i)
        getvecmapLinkerorMotor<LinkerorMotor>()[i].clear();

//    cout<<"size of dense hash map "<<N<<endl;

    //Check if there are entries
/*
    for(auto key: googlepossiblel){
        cout<<key.first<<" "<<key.second.size()<<endl;
    }
*/

}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::threadwritepairsself(uint startID, uint endID,
     dist::dOut<D, SELF> &bspairsoutS, short idvec[2],
     gdmap tempmap[], short threadID) {

    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> rtemp;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);

    auto cylcmp1 = _compartment->Cyldcindexvec;
    short idx = idvec[0];
//    short idx2 = idvec[1];
//    short Nbounds = _rMaxsqvec[idx].size();
    short dim = 0;
    for(uint pid = startID; pid < endID; pid++) {

        uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;

        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        int cIndex1 = cylcmp1[cmpcylidx1];
        int cIndex2 = cylcmp1[cmpcylidx2];

        auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

        auto fpairs = _filamentIDvec[idx].data();
        //Check if other conditions are satisfied
        if(cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
           .filamentposition - cylinder2.filamentposition) <=2)
            continue;
        if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] || cylinder1
             .type == fpairs[1] && cylinder2.type == fpairs[0])) {

            //for(uint idx = 0; idx < N; idx++){
            if (cIndex1 > cIndex2) {

                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                temp[t1].push_back(t2);
                rtemp[t2].push_back(t1);
            } else {
                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                temp[t2].push_back(t1);
                rtemp[t1].push_back(t2);
            }
        }
    }
    tempmap[2 * threadID] = temp;
    tempmap[2 * threadID + 1] = rtemp;
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::threadwritepairs(Compartment* ncmp,
        dist::dOut<D, SELF> &bspairsoutS,
        gdmap tempmap[], gdmap tempNmap[],
        uint startID, uint endID, short idvec[2], short threadID) {

    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> rtemp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> ntemp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> nrtemp;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);
    ntemp.set_empty_key(4294967295);
    ntemp.set_deleted_key(4294967294);
    nrtemp.set_empty_key(4294967295);
    nrtemp.set_deleted_key(4294967294);

    auto cylcmp1 = _compartment->Cyldcindexvec;
    auto cylcmp2 = ncmp->Cyldcindexvec;

    short idx = idvec[0];
//    short idx2 = idvec[1];
//    short Nbounds = _rMaxsqvec[idx].size();
    short dim = 0;
    for(uint pid = startID; pid < endID; pid++) {

        uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;

        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        int cIndex1 = cylcmp1[cmpcylidx1];
        int cIndex2 = cylcmp2[cmpcylidx2];

        auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

        auto fpairs = _filamentIDvec[idx].data();

        //Check if other conditions are satisfied
        if(cylinder1.filamentID == cylinder2.filamentID &&
            abs(cylinder1.filamentposition - cylinder2.filamentposition) <=2)
            continue;

        if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] ||
            cylinder1.type == fpairs[1] && cylinder2.type == fpairs[0])) {

            //for(uint idx = 0; idx < N; idx++){
            if (cIndex1 > cIndex2) {

                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                temp[t1].push_back(t2);
                rtemp[t2].push_back(t1);
            } else {
                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                ntemp[t2].push_back(t1);
                nrtemp[t1].push_back(t2);
            }
        }
    }

    tempmap[2 * threadID] = temp;
    tempmap[2 * threadID + 1] = rtemp;
    tempNmap[2 * threadID] = ntemp;
    tempNmap[2 * threadID + 1] = nrtemp;
}

//SECOND VERSION V2
template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::threadwritepairsselfV2(uint startID, uint endID, dist::dOut<D,SELF>& bspairsoutS,
                                                        short idvec[2],
                                                        dense_hash_map<uint32_t , uint, hash<uint32_t>> tempmap[],
                                                        vector<vector<uint32_t>>
                                                                valuematrixvec[],
                                                        short threadID) {

    dense_hash_map<uint32_t, uint32_t, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, uint32_t, hash<uint32_t>> rtemp;
    vector<vector<uint>> tempmatrix;
    vector<vector<uint>> rtempmatrix;

    uint rowID = 0;
    uint rowIDr = 0;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);

    auto cylcmp1 = _compartment->Cyldcindexvec;
    short idx = idvec[0];
//    short idx2 = idvec[1];
//    short Nbounds = _rMaxsqvec[idx].size();
    short dim = 0;
    for(uint pid = startID; pid < endID; pid++) {

        uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;

        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        int cIndex1 = cylcmp1[cmpcylidx1];
        int cIndex2 = cylcmp1[cmpcylidx2];

        auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

        auto fpairs = _filamentIDvec[idx].data();
        //Check if other conditions are satisfied
        if(cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1
           .filamentposition - cylinder2.filamentposition) <=2)
            continue;
        if ((cylinder1.type == fpairs[0] && cylinder2.type == fpairs[1] || cylinder1
             .type == fpairs[1] && cylinder2.type == fpairs[0])) {

            //for(uint idx = 0; idx < N; idx++){
            if (cIndex1 > cIndex2) {

                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;

                if(temp.find(t1) == temp.end()){
                    temp[t1] = rowID;
                    rowID++;
                }

                if(rtemp.find(t2) == rtemp.end()){
                    rtemp[t2] = rowIDr;
                    rowIDr++;
                }

                if(tempmatrix.size()==temp[t1])
                    tempmatrix.push_back({t2});
                else
                    tempmatrix.at(temp[t1]).push_back(t2);

                if(rtempmatrix.size() == rtemp[t2])
                    rtempmatrix.push_back({t1});
                else
                    rtempmatrix.at(rtemp[t2]).push_back(t1);

            } else {

                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;

                if(temp.find(t2) == temp.end()){
                    temp[t2] = rowID;
                    rowID++;
                }

                if(rtemp.find(t1) == rtemp.end()){
                    rtemp[t1] = rowIDr;
                    rowIDr++;
                }
                if(tempmatrix.size()==temp[t2])
                    tempmatrix.push_back({t1});
                else
                    tempmatrix.at(temp[t2]).push_back(t1);

                if(rtempmatrix.size() == rtemp[t1])
                    rtempmatrix.push_back({t2});
                else
                    rtempmatrix.at(rtemp[t1]).push_back(t2);
            }
        }
    }
    tempmap[2 * threadID] = temp;
    tempmap[2 * threadID + 1] = rtemp;
    valuematrixvec[2 * threadID] = tempmatrix;
    valuematrixvec[2 * threadID + 1] = rtempmatrix;
}

//THIRD VERSION V3
template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::threadwritepairsselfV3(uint startID, uint endID,
                                                      dist::dOut<D, SELF> &bspairsoutS, short idvec[2],
                                                      short threadID) {

    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> rtemp;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);

    auto cylcmp1 = _compartment->Cyldcindexvec;
    short idx = idvec[0];
//    short idx2 = idvec[1];
//    short Nbounds = _rMaxsqvec[idx].size();
    short dim = 0;
    for(uint pid = startID; pid < endID; pid++) {

        uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;

        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        int cIndex1 = cylcmp1[cmpcylidx1];
        int cIndex2 = cylcmp1[cmpcylidx2];

/*        auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];*/

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

        auto fpairs = _filamentIDvec[idx].data();
        //Check if other conditions are satisfied
/*        if(cylinder1.filamentID == cylinder2.filamentID && abs(cylinder1.filamentposition -
        cylinder2.filamentposition) <=2)
            continue;*/
        {

                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                temp[t1].push_back(t2);
                rtemp[t2].push_back(t1);
        }
    }

    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID] = temp;
    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID + 1] = rtemp;
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::threadwritepairsV3(Compartment* ncmp,
                                                  dist::dOut<D, SELF> &bspairsoutS,
                                                  uint startID, uint endID, short idvec[2],
                                                  short threadID) {

    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> rtemp;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);

    auto cylcmp1 = _compartment->Cyldcindexvec;
    auto cylcmp2 = ncmp->Cyldcindexvec;

    short idx = idvec[0];
//    short idx2 = idvec[1];
//    short Nbounds = _rMaxsqvec[idx].size();
    short dim = 0;
    for(uint pid = startID; pid < endID; pid++) {

        uint16_t site1 = bspairsoutS.dout[2 * dim][pid];
        uint16_t site2 = bspairsoutS.dout[2 * dim + 1][pid];

        uint16_t cmpcylidx1 = site1;
        uint16_t cmpcylidx2 = site2;

        cmpcylidx1 = (cmpcylidx1 >> 4);
        cmpcylidx2 = (cmpcylidx2 >> 4);

        int cIndex1 = cylcmp1[cmpcylidx1];
        int cIndex2 = cylcmp2[cmpcylidx2];

/*        auto cylinder1 = CUDAcommon::serlvars.cylindervec[cIndex1];
        auto cylinder2 = CUDAcommon::serlvars.cylindervec[cIndex2];*/

        short bsite1 = mask & site1;
        short bsite2 = mask & site2;

        auto fpairs = _filamentIDvec[idx].data();

        //Check if other conditions are satisfied
/*        if(cylinder1.filamentID == cylinder2.filamentID &&
           abs(cylinder1.filamentposition - cylinder2.filamentposition) <=2)
            continue;*/
        {
                uint32_t t1 = cIndex1<<4|bsite1;
                uint32_t t2 = cIndex2<<4|bsite2;
                temp[t1].push_back(t2);
                rtemp[t2].push_back(t1);
        }
    }

    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID] = temp;
    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID + 1] = rtemp;
}

template<uint D, bool LinkerorMotor>
void HybridBindingSearchManager::threadedsingleparse(uint first, uint last, short threadID){

    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> temp;
    dense_hash_map<uint32_t, vector<uint32_t>, hash<uint32_t>> rtemp;

    temp.set_empty_key(4294967295);
    temp.set_deleted_key(4294967294);
    rtemp.set_empty_key(4294967295);
    rtemp.set_deleted_key(4294967294);

    auto LorMpairs = getpairsLinkerorMotor<LinkerorMotor>();

//    auto cylcmp1 = _compartment->Cyldcindexvec;

    for(uint pid = first; pid < last; pid++) {


        //Check if other conditions are satisfied
        if(getfilID_fpospairs<LinkerorMotor>()[pid]) {
            uint32_t t1 = LorMpairs[2 * pid];
            uint32_t t2 = LorMpairs[2 * pid + 1];
            temp[t1].push_back(t2);
            rtemp[t2].push_back(t1);
        }
/*        else{
            cout<<"OOPS! Mistake "<<endl;
        }*/
    }


    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID] = temp;
    getvecmapLinkerorMotor<LinkerorMotor>()[2 * threadID + 1] = rtemp;

}

template<bool LinkerorMotor>
void HybridBindingSearchManager::mergemaps(short threadID, short offset){

    short parentID = threadID;
    short partnerID = (threadID + offset);
    if (getvecmapLinkerorMotor<LinkerorMotor>()[partnerID].size() ) {
        for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[partnerID]) {

            vector<uint32_t> sitestoappend = iter.second;
            auto mapvalue = getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first];
//            uint bfz = mapvalue.size();
            getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].reserve
                    (mapvalue.size() + sitestoappend.size());
            getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].insert
                    (getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].end(),
                     sitestoappend.begin(), sitestoappend.end());
/*            if (mapvalue.size() == 0) {
                sitestoappend.reserve(sitestoappend.size()+1024);
                getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].insert
                (getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].end(),
                        sitestoappend.begin(), sitestoappend.end());
            } else {
                getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].reserve
                (mapvalue.size() + sitestoappend.size());
                getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].insert
                (getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].end(),
                        sitestoappend.begin(), sitestoappend.end());
            }*/
/*            uint afz = getvecmapLinkerorMotor<LinkerorMotor>()[parentID][iter.first].size();
            if(sitestoappend.size() == 0 || afz ==0)
                cout<<"Key "<<iter.first<<" Before "<<bfz<<" After "<<afz<<" additional size "
                <<sitestoappend.size()<<endl;*/

            /*auto mapiter = getvecmapLinkerorMotor<LinkerorMotor>()[parentID].find(iter.first);
            if (mapiter == getvecmapLinkerorMotor<LinkerorMotor>()[parentID].end()) {
                sitestoappend.reserve(sitestoappend.size()+1024);
                getvecmapLinkerorMotor<LinkerorMotor>()[parentID].insert(make_pair(iter.first, sitestoappend));
            } else {
                mapiter->second.resize(mapiter->second.size() + sitestoappend.size());
                mapiter->second.insert(mapiter->second.end(), sitestoappend.begin(),
                                       sitestoappend.end());
            }*/
        }
    }
    //reverse map
    /*if (getvecmapLinkerorMotor<LinkerorMotor>()[partnerID + 1].size()) {
        for (const auto &iter : getvecmapLinkerorMotor<LinkerorMotor>()[partnerID + 1]) {

            vector<uint32_t> sitestoappend = iter.second;
            auto mapiter = getvecmapLinkerorMotor<LinkerorMotor>()[parentID + 1].find(iter.first);
            if (mapiter == getvecmapLinkerorMotor<LinkerorMotor>()[parentID + 1].end()) {
                sitestoappend.reserve(sitestoappend.size()+1024);
                getvecmapLinkerorMotor<LinkerorMotor>()[parentID + 1].insert(make_pair(iter.first, sitestoappend));
            } else {
                mapiter->second.resize(mapiter->second.size() + sitestoappend.size());
                mapiter->second.insert(mapiter->second.end(), sitestoappend.begin(),
                                       sitestoappend.end());
            }
        }
    }*/
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
        for(short idx2 = 0; idx2<temprMaxsq.size();idx2++) {
                        short temp = 0;
#ifndef SIMDBINDINGSEARCH
            temp = _HneighborList->setneighborsearchparameters(ftypepair[0], ftypepair[1],
                          false, true, localmaxcylsize + sqrt(_rMaxsqvec[idx][idx2]),
                          max( sqrt(_rMinsqvec[idx][idx2]) - localmaxcylsize,float(0.0)));
#endif
            HNLIDvec[idx]=(temp);
            short idvec[2] = {idx, idx2};
            fManagervec[idx][idx2]->setHNLID(HNLIDvec[idx], idvec);
        }
    }
    //Set google vars
/*    if(!googlevar){
        short idx = 0;
        for(short idx2 = 0;idx2<2;idx2++) {
            googlepossible[idx][idx2].set_empty_key(4294967295);
            googlepossible[idx][idx2].set_deleted_key(4294967294);
            googlereversepossible[idx][idx2].set_empty_key(4294967295);
            googlereversepossible[idx][idx2].set_deleted_key(4294967294);
        }
        googlevar = true;
    }*/

}

vector<tuple<CCylinder*, short>>
HybridBindingSearchManager::chooseBindingSitesstencil(short idvec[2]){
#ifdef SIMD_BINDINGSEARCH
    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fpairs = _filamentIDvec[idx].data();
    int pbsSize = Nbindingpairs[idx][idx2];
    assert((pbsSize!= 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");
    uint randomIndex = Rand::randInteger(1, pbsSize);

    auto it = googlepossible[idx][idx2].begin();
    uint cumulativesum = 0;
    uint prevstepsum = 0;

    while(it != googlepossible[idx][idx2].end() && cumulativesum < randomIndex){
        prevstepsum = cumulativesum;
        cumulativesum += it->second.size();
        if(cumulativesum < randomIndex)
            it++;
    }
    Nbindingpairs[idx][idx2]--;

    uint position  = randomIndex - prevstepsum -1;

    uint32_t site1 = it->first;
    uint32_t site2 = it->second[position];

    uint32_t CylID1 = site1>>4;
    uint32_t CylID2 = site2>>4;

    short bsitepos1 = mask & site1;
    short bsitepos2 = mask & site2;

/*    if(bsitepos1>3 || bsitepos2>3){
        cout<<"check possible "<<endl;
        cout<<"Site1 "<<bitset<32>(site1)<<" CylID "<<bitset<32>(CylID1)<<endl;
        cout<<"Site2 "<<bitset<32>(site2)<<" CylID "<<bitset<32>(CylID2)<<endl;
    }*/

    bool foundcyl1 = false;
    bool foundcyl2 = false;

    CCylinder* ccyl1;
    CCylinder* ccyl2;

    auto Cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    cout<<"CHECK "<<CylID1<<" "<<CylID2<<" "<<position<<" "<<it->second.size()<<endl;
    uint i = 0;
    while(!(foundcyl1 && foundcyl2)){
        if(i==Cylinder::getCylinders().size()){
            cout<<"Error! Couldn't find CylinderID"<<endl;
        }
        uint ID = Cylinderpointervec[i]->getID();
        if(ID == CylID1){
            ccyl1 = Cylinderpointervec[i]->getCCylinder();
            fpairs[0] = Cylinderpointervec[i]->getType();
            foundcyl1 = true;
//            cout<<"ID1 "<<ID<<" "<<ccyl1->getCylinder()->getID()<<endl;
        }
        else if(ID == CylID2){
            ccyl2 = Cylinderpointervec[i]->getCCylinder();
            fpairs[1] =  Cylinderpointervec[i]->getType();
            foundcyl2 = true;
//            cout<<"ID2 "<<ID<<" "<<ccyl2->getCylinder()->getID()<<endl;
        }
        i++;
    }
//    cout<<"chosen CylID "<<CylID1<<" "<<CylID2<<" bs "<<bsitepos1<<" "<<bsitepos2<<" fID "
//    <<fpairs[0]<<" "<<fpairs[1]<<endl;
    short bindingSite1 = SysParams::Chemistry().bindingSites[fpairs[0]][bsitepos1];
    short bindingSite2 = SysParams::Chemistry().bindingSites[fpairs[1]][bsitepos2];

    tuple<CCylinder*, short> t1 = make_tuple(ccyl1, bindingSite1);
    tuple<CCylinder*, short> t2 = make_tuple(ccyl2, bindingSite2);

    return vector<tuple<CCylinder*, short>>{t1, t2};
#else
    short idx = idvec[0];
    short idx2 = idvec[1];
    int pbsSize = _possibleBindingsstencilvec[idx][idx2].size() ;
    cout<<idvec[0]<<" "<<idvec[1]<<" "<<_possibleBindingsstencilvec[idx][idx2].size()<<endl;
    assert((pbsSize!= 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");
    int randomIndex = Rand::randInteger(0, pbsSize - 1);
    auto it = _possibleBindingsstencilvec[idx][idx2].begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};

#endif

}

HybridCylinderCylinderNL* HybridBindingSearchManager::_HneighborList;
bool initialized = false;
//D = 1
dist::dOut<1U,true> HybridBindingSearchManager::bspairs1self;
dist::dOut<1U,false> HybridBindingSearchManager::bspairs1;
dist::dOut<1U,false> HybridBindingSearchManager::bspairslinker;
dist::dOut<1U,true> HybridBindingSearchManager::bspairslinkerself;
dist::dOut<1U,false> HybridBindingSearchManager::bspairsmotor;
dist::dOut<1U,true> HybridBindingSearchManager::bspairsmotorself;
//D = 2
dist::dOut<2U,true> HybridBindingSearchManager::bspairs2self;
dist::dOut<2U,false> HybridBindingSearchManager::bspairs2;
/*dist::dOut<1U,false> HybridBindingSearchManager::bspairs2_D1;
dist::dOut<1U,false> HybridBindingSearchManager::bspairs2_D2;*/


short HybridBindingSearchManager::Totallinkermotor = 0;
double HybridBindingSearchManager::SIMDtime = 0.0;
double HybridBindingSearchManager::HYBDtime = 0.0;
double HybridBindingSearchManager::findtime = 0.0;
double HybridBindingSearchManager::appendtime = 0.0;
double HybridBindingSearchManager::findtimeV2 = 0.0;
double HybridBindingSearchManager::SIMDparse1 = 0.0;
double HybridBindingSearchManager::SIMDparse2 = 0.0;
double HybridBindingSearchManager::SIMDparse3 = 0.0;
double HybridBindingSearchManager::SIMDcountbs = 0.0;


/*///Template specializations
template void HybridBindingSearchManager::calculatebspairsself<1,true>();
template void HybridBindingSearchManager::calculatebspairsself<1,false>();
template void HybridBindingSearchManager::calculatebspairsself<2,true>();
template void HybridBindingSearchManager::calculatebspairsself<2,false>();
template void HybridBindingSearchManager::calculatebspairsself<3,true>();
template void HybridBindingSearchManager::calculatebspairsself<3,false>();
template void HybridBindingSearchManager::calculatebspairsself<4,true>();
template void HybridBindingSearchManager::calculatebspairsself<4,false>();*/
#endif

