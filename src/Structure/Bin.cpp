
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
#include "Bin.h"
#include "Cylinder.h"
void Bin::updatecindices(){
    cindicesvector.clear();
    cindicesvector.reserve(_cylinders.size());
    for(auto &c:_cylinders)
        cindicesvector.push_back(c->_dcIndex);
}

#ifdef SIMDBINDINGSEARCH2
void Bin::createSIMDcoordinates(){
    int N = _cylinders.size() * SysParams::Chemistry().maxbindingsitespercylinder;
//    cout<<"BinID "<<_ID<<" #cyls "<<N<<endl;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(
                N), bindsitecoordinatesZ(N);
        cindex_bs.resize(N);
        //cID_bs.resize(N);
        Cyldcindexvec.resize(_cylinders.size());
        CylcIDvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        for (auto cyl:_cylinders) {
            int cindex = cyl->_dcIndex;
            auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
            if (checkftype)
                short _filamentType = cylinderstruct.type;
            auto x1 = cyl->getFirstBead()->coordinate;
            auto x2 = cyl->getSecondBead()->coordinate;
            uint16_t shiftedindex = (i << 4);
            //uint32_t shiftedCID = cyl->getID()<<4;
            Cyldcindexvec[i] = cyl->_dcIndex;
            CylcIDvec[i] = cyl->getID();
            i++;
            uint16_t j = 0;
            for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                 it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {

                auto mp = (float) *it / SysParams::Geometry().cylinderNumMon[_filamentType];
                auto coord = midPointCoordinate(x1, x2, mp);
                bindsitecoordinatesX[k] = coord[0];
                bindsitecoordinatesY[k] = coord[1];
                bindsitecoordinatesZ[k] = coord[2];
                //last 4 bits are binding site while first 12 bits are cylinder index.
                cindex_bs[k] = shiftedindex | j;
                //cID_bs[k] = shiftedCID |j;
                k++;
                j++;

            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
        bscoords.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                             bindsitecoordinatesZ, cindex_bs);
    }

}

void Bin::createSIMDcoordinates4linkersearch(bool isvectorizedgather){
//setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 1;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    int N = _cylinders.size() * maxnbs;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(N),
                bindsitecoordinatesZ(N);
        vector<uint16_t> cindex_bs(N);
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        if(isvectorizedgather) {
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (boundstate[bstatepos][maxnbs * cindex + j]) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
        else{
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if(cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN
                            () == 1.0) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
        bscoordslinker.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                                   bindsitecoordinatesZ,
                                   cindex_bs);
    }
}

void Bin::createSIMDcoordinates4motorsearch(bool isvectorizedgather){
//setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 2;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    int N = _cylinders.size() * maxnbs;
    if(N) {
        vector<double> bindsitecoordinatesX(N), bindsitecoordinatesY(N),
                bindsitecoordinatesZ(N);
        vector<uint16_t> cindex_bs(N);
        Cyldcindexvec.resize(_cylinders.size());

        short _filamentType = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;
        uint16_t k = 0;
        if(isvectorizedgather) {
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (boundstate[bstatepos][maxnbs * cindex + j]) {
                        int A = boundstate[bstatepos][maxnbs * cindex + j] -
                                cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                                        SysParams::Chemistry().motorBoundIndex[_filamentType])->getN
                                        ();
                        if(abs(A) != 0){
                            cout<<"OOPS Problematic! Occupies sites are passed to "
                                  "SIMD"<<endl;
                            cout<<"Cylinder "<<cyl->getID()<<" "<<*it<<endl;
                        }
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
        else{
            for (auto cyl:_cylinders) {
                int cindex = cyl->_dcIndex;
                auto cylinderstruct = CUDAcommon::serlvars.cylindervec[cindex];
                if (checkftype)
                    short _filamentType = cylinderstruct.type;
                auto x1 = cyl->getFirstBead()->coordinate;
                auto x2 = cyl->getSecondBead()->coordinate;
                uint16_t shiftedindex = (i << 4);
                Cyldcindexvec[i] = cindex;
                i++;
                uint16_t j = 0;
                for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                     it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                    if (cyl->getCCylinder()->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN
                            () == 1.0) {
                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto coord = midPointCoordinate(x1, x2, mp);
                        bindsitecoordinatesX[k] = coord[0];
                        bindsitecoordinatesY[k] = coord[1];
                        bindsitecoordinatesZ[k] = coord[2];
                        //last 4 bits are binding site while first 12 bits are cylinder index.
                        cindex_bs[k] = shiftedindex | j;
                        k++;
                    }
                    j++;
                }
            }
        }
//    assert(k<65536);
        //Create input vector for SIMD calculations
        bscoordsmotor.init_coords(bindsitecoordinatesX, bindsitecoordinatesY,
                                  bindsitecoordinatesZ,
                                  cindex_bs);
    }
}
#endif